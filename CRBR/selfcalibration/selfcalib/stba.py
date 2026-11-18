"""Spatiotemporal self-calibrating bundle adjustment (STBA) from manual labels.

Model
-----
Camera c's label row r (0-based) is captured at absolute time

    t = (r + 1 - off_c) - d_c  [+ rho * (v_raw / H - 0.5)]

where off_c is the integer offset already in offset.xlsx, d_c a continuous
correction in frames (the same convention as sync.py: d_c is ADDED to the
offset), and rho an optional global rolling-shutter readout (frames per image
height).  Each landmark j is a free 3-D trajectory X_j on the integer frame grid,
evaluated at fractional t with Catmull-Rom interpolation.  The label is
predicted with the full fisheye camera model (raw, distorted MATLAB pixels):

    u = pi_c( X_j(t) ),   camera c: R_c = exp([w_c]x) R_c0,  C_c = C_c0 + dC_c,
                           K_c = K_c0 with (fx, fy) * exp(s_c)  (optional)

Cost = sum Huber(reprojection / sigma) + sum Huber(2nd-difference / sigma_acc)
       + Gaussian priors on w, dC, s, d  (the priors also fix the gauge).

Free-point structure (no rigid template) keeps camera/sync estimation
independent of the anatomical model used later by the CRBR reconstruction.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix

from .camera import Camera, kb_distort, rodrigues

CALIB_POINTS = ["l_eye", "r_eye", "nose", "l_ear", "r_ear", "spine", "neck_base", "l_shoulder", "r_shoulder"]


def catmull_rom_weights(u):
    """Weights for points (i-1, i, i+1, i+2) at local parameter u in [0,1)."""
    u2, u3 = u * u, u * u * u
    return 0.5 * np.stack([-u + 2 * u2 - u3,
                           2 - 5 * u2 + 3 * u3,
                           u + 4 * u2 - 3 * u3,
                           -u2 + u3], axis=-1)


class STBA:
    def __init__(self, ds, points=None, cams=None, d0=None, refine_calib=True, refine_sync=True,
                 refine_focal=False, rolling_shutter=False, ref_cam=None,
                 sigma_rot_deg=1.0, sigma_pos_m=0.10, sigma_focal=0.01, sigma_sync=2.0, sync_bound=1.5,
                 acc_ms2=100.0, fps=90.0, huber_reproj=2.5, huber_acc=3.0, obs_filter=None,
                 drop_cams_obs=(), margin=3, sync_cams=None, gauge="similarity", calib_min_obs=0,
                 freeze_calib_cams=()):
        self.ds = ds
        self.points = [p for p in (CALIB_POINTS if points is None else points) if p in ds.points]
        self.cams0 = [c.copy() for c in (ds.cams if cams is None else cams)]
        C = ds.ncams
        self.C = C
        self.d0 = np.zeros(C) if d0 is None else np.asarray(d0, float).copy()
        self.refine_calib, self.refine_sync = refine_calib, refine_sync
        self.refine_focal, self.rolling_shutter = refine_focal, rolling_shutter
        self.sync_cams = None if sync_cams is None else set(int(k) for k in sync_cams)
        # sync prior is anchored on the ORIGINAL offsets (d = 0), not on the
        # initialiser: a camera whose data cannot inform its sync falls back to
        # the manual offset instead of drifting to an arbitrary start value.
        self.d_prior = np.zeros(ds.ncams)
        self.gauge, self.calib_min_obs = gauge, calib_min_obs
        self.freeze_calib_cams = set(int(k) for k in freeze_calib_cams)
        self.sig_rot = np.deg2rad(sigma_rot_deg)
        self.sig_pos, self.sig_focal, self.sig_sync = sigma_pos_m, sigma_focal, sigma_sync
        self.sync_bound = sync_bound
        self.sig_acc = acc_ms2 / fps ** 2
        self.hub_r, self.hub_a = huber_reproj, huber_acc
        self.H = ds.H

        # ---- observations ----
        raws, _, sigs = ds.camera_tracks(self.cams0)
        rows_c, cam_c, pt_c, u_c, s_c = [], [], [], [], []
        for k in range(C):
            if k in drop_cams_obs:
                continue
            for jj, p in enumerate(self.points):
                j = ds.pidx(p)
                r = np.where(np.all(np.isfinite(raws[k][:, j]), -1))[0]
                if obs_filter is not None:
                    tb = r + 1 - ds.offsets[k]
                    r = r[obs_filter(k, p, tb)]
                rows_c.append(r); cam_c.append(np.full(len(r), k)); pt_c.append(np.full(len(r), jj))
                u_c.append(raws[k][r, j]); s_c.append(sigs[k][r, j])
        self.o_row = np.concatenate(rows_c).astype(float)
        self.o_cam = np.concatenate(cam_c).astype(int)
        self.o_pt = np.concatenate(pt_c).astype(int)
        self.o_u = np.concatenate(u_c)
        self.o_sig = np.concatenate(s_c)
        self.o_tb = self.o_row + 1 - ds.offsets[self.o_cam]
        self.o_rs = self.o_u[:, 1] / self.H - 0.5
        self.nobs = len(self.o_row)

        # cameras that actually have observations
        self.cam_has = np.array([np.any(self.o_cam == k) for k in range(C)])
        if ref_cam is None:
            cnt = np.bincount(self.o_cam, minlength=C)
            ref_cam = int(np.argmax(cnt))
        self.ref = ref_cam

        # ---- reachable sync range per camera (matches bounds()) ----
        rs = 1.0 if rolling_shutter else 0.0
        if refine_sync:
            self.d_lo = np.minimum(self.d0, self.d_prior) - sync_bound - rs
            self.d_hi = np.maximum(self.d0, self.d_prior) + sync_bound + rs
        else:
            self.d_lo = self.d0 - rs
            self.d_hi = self.d0 + rs
        # ---- trajectory grids per point ----
        self.grid = []
        for jj in range(len(self.points)):
            m = self.o_pt == jj
            if not m.any():
                self.grid.append((0, -1)); continue
            kk = self.o_cam[m]
            g0 = int(np.floor(np.min(self.o_tb[m] - self.d_hi[kk]))) - margin
            g1 = int(np.ceil(np.max(self.o_tb[m] - self.d_lo[kk]))) + margin
            self.grid.append((g0, g1))
        self._layout()
        self.x = self._init_x()
        self.w_obs = np.ones(self.nobs)
        self.w_acc = None

    # ------------------------------------------------------------------
    def _layout(self):
        idx = 0
        self.cam_slices = {}
        for k in range(self.C):
            if not self.cam_has[k]:
                continue
            sl = {}
            if (self.refine_calib and np.sum(self.o_cam == k) >= self.calib_min_obs
                    and k not in self.freeze_calib_cams):
                sl["w"] = np.arange(idx, idx + 3); idx += 3
                sl["dC"] = np.arange(idx, idx + 3); idx += 3
            if self.refine_focal:
                sl["s"] = np.arange(idx, idx + 1); idx += 1
            if self.refine_sync and k != self.ref and (self.sync_cams is None or k in self.sync_cams):
                sl["d"] = np.arange(idx, idx + 1); idx += 1
            self.cam_slices[k] = sl
        self.rho_idx = None
        if self.rolling_shutter:
            self.rho_idx = idx; idx += 1
        self.ncam_params = idx
        self._gauge_basis()
        self.traj_off = []
        for (g0, g1) in self.grid:
            n = max(0, g1 - g0 + 1)
            self.traj_off.append(idx)
            idx += 3 * n
        self.nparams = idx

    def _init_x(self):
        from .triangulate import triangulate
        x = np.zeros(self.nparams)
        for k, sl in self.cam_slices.items():
            if "d" in sl:
                x[sl["d"]] = self.d0[k]
        ds = self.ds
        # structure init: triangulate at integer frames from labels nearest in time
        for jj, (g0, g1) in enumerate(self.grid):
            if g1 < g0:
                continue
            n = g1 - g0 + 1
            m = self.o_pt == jj
            tt = self.o_tb[m] - self.d0[self.o_cam[m]]
            P = np.full((n, 3), np.nan)
            for gi in range(n):
                f = g0 + gi
                sel = np.where(m)[0][np.abs(tt - f) < 0.5 + 1e-9]
                ks = self.o_cam[sel]
                if len(np.unique(ks)) < 2:
                    continue
                _, first = np.unique(ks, return_index=True)
                sel = sel[first]
                cams = [self.cams0[k] for k in self.o_cam[sel]]
                uv = np.stack([cams[i].undistort(self.o_u[s]) for i, s in enumerate(sel)])
                X, _, dep = triangulate(cams, uv, w=1.0 / self.o_sig[sel])
                if np.all(np.isfinite(X)) and np.all(dep > 0):
                    P[gi] = X
            ok = np.all(np.isfinite(P), 1)
            if not ok.any():
                # fall back: centre of all cameras' view (rough)
                P[:] = np.mean([c.center for c in self.cams0], 0)
            else:
                fi = np.arange(n)
                for a in range(3):
                    P[:, a] = np.interp(fi, fi[ok], P[ok, a])
            x[self.traj_off[jj]: self.traj_off[jj] + 3 * n] = P.ravel()
        return x

    # ------------------------------------------------------------------
    def cameras(self, x=None):
        x = self.x if x is None else x
        out = []
        for k in range(self.C):
            c = self.cams0[k].copy()
            sl = self.cam_slices.get(k, {})
            if "w" in sl:
                R = rodrigues(x[sl["w"]]) @ self.cams0[k].R
                Cc = self.cams0[k].center + x[sl["dC"]]
                c.R = R
                c.t = -R @ Cc
            if "s" in sl:
                e = np.exp(x[sl["s"]][0])
                c.K = c.K.copy(); c.K[0, 0] *= e; c.K[1, 1] *= e
            out.append(c)
        return out

    def offsets_corr(self, x=None):
        x = self.x if x is None else x
        d = self.d0.copy()
        for k, sl in self.cam_slices.items():
            if "d" in sl:
                d[k] = x[sl["d"]][0]
        return d

    def rho(self, x=None):
        x = self.x if x is None else x
        return x[self.rho_idx] if self.rho_idx is not None else 0.0

    def obs_times(self, x=None):
        d = self.offsets_corr(x)
        return self.o_tb - d[self.o_cam] + self.rho(x) * self.o_rs

    def traj(self, jj, x=None):
        x = self.x if x is None else x
        g0, g1 = self.grid[jj]
        n = g1 - g0 + 1
        return x[self.traj_off[jj]: self.traj_off[jj] + 3 * n].reshape(n, 3)

    def eval_points(self, x, t, jj_arr):
        """Interpolated 3-D positions for observation times t of points jj_arr."""
        X = np.zeros((len(t), 3))
        for jj in np.unique(jj_arr):
            m = jj_arr == jj
            g0, g1 = self.grid[jj]
            P = self.traj(jj, x)
            n = len(P)
            tt = t[m] - g0
            i = np.clip(np.floor(tt).astype(int), 1, n - 3)
            u = tt - i
            W = catmull_rom_weights(u)
            X[m] = (W[:, 0:1] * P[i - 1] + W[:, 1:2] * P[i] + W[:, 2:3] * P[i + 1] + W[:, 3:4] * P[i + 2])
        return X

    def project_obs(self, x):
        cams = self.cameras(x)
        t = self.obs_times(x)
        X = self.eval_points(x, t, self.o_pt)
        uv = np.zeros((self.nobs, 2))
        for k in range(self.C):
            m = self.o_cam == k
            if m.any():
                uv[m] = cams[k].project(X[m], distort=True)
        return uv

    def reproj_residuals(self, x=None):
        x = self.x if x is None else x
        return (self.project_obs(x) - self.o_u) / self.o_sig[:, None]

    def acc_residuals(self, x):
        out = []
        for jj in range(len(self.points)):
            g0, g1 = self.grid[jj]
            if g1 - g0 < 2:
                continue
            P = self.traj(jj, x)
            out.append((P[2:] - 2 * P[1:-1] + P[:-2]) / self.sig_acc)
        return np.concatenate(out) if out else np.zeros((0, 3))

    def prior_residuals(self, x):
        r = []
        for k, sl in self.cam_slices.items():
            if "w" in sl:
                r.append(x[sl["w"]] / self.sig_rot)
                r.append(x[sl["dC"]] / self.sig_pos)
            if "s" in sl:
                r.append(x[sl["s"]] / self.sig_focal)
            if "d" in sl:
                r.append((x[sl["d"]] - self.d_prior[k]) / self.sig_sync)
        if self.rho_idx is not None:
            r.append(np.array([x[self.rho_idx] / 2.0]))
        if self.gbasis is not None:
            r.append(self.gbasis @ x[self.gcols] * 100.0)
        return np.concatenate(r) if r else np.zeros(0)

    def _gauge_basis(self):
        """Rows = the 7 infinitesimal similarity transforms of the whole rig,
        expressed in the (prior-whitened) camera-correction coordinates.  The
        residual G x = 0 removes the common-mode (gauge) component of the
        corrections, so a weakly observed camera that stays near its original
        pose remains consistent with the refined ones."""
        self.gbasis, self.gcols = None, None
        cal = [k for k, sl in self.cam_slices.items() if "w" in sl]
        if self.gauge != "similarity" or len(cal) < 2:
            return
        cols = np.concatenate([np.r_[self.cam_slices[k]["w"], self.cam_slices[k]["dC"]] for k in cal])
        Cs = np.array([self.cams0[k].center for k in cal])
        Cm = Cs.mean(0)
        B = []
        E = np.eye(3)
        for a in range(3):  # rotations
            v = []
            for i, k in enumerate(cal):
                v.append(np.r_[-self.cams0[k].R @ E[a] / self.sig_rot ** 2, np.cross(E[a], Cs[i] - Cm) / self.sig_pos ** 2])
            B.append(np.concatenate(v))
        for a in range(3):  # translations
            B.append(np.concatenate([np.r_[np.zeros(3), E[a] / self.sig_pos ** 2] for k in cal]))
        B.append(np.concatenate([np.r_[np.zeros(3), (Cs[i] - Cm) / self.sig_pos ** 2] for i, k in enumerate(cal)]))
        B = np.array(B)
        # normalise each constraint to unit norm in whitened units
        W = np.concatenate([np.r_[np.full(3, self.sig_rot), np.full(3, self.sig_pos)] for k in cal])
        B = B / np.linalg.norm(B * W[None, :], axis=1, keepdims=True)
        self.gbasis, self.gcols = B, cols

    def fun(self, x):
        rr = self.reproj_residuals(x) * np.sqrt(self.w_obs)[:, None]
        ra = self.acc_residuals(x)
        if self.w_acc is not None:
            ra = ra * np.sqrt(self.w_acc)[:, None]
        return np.concatenate([rr.ravel(), ra.ravel(), self.prior_residuals(x)])

    # ------------------------------------------------------------------
    def sparsity(self, free_cam=True):
        nr = 2 * self.nobs + 3 * sum(max(0, g1 - g0 - 1) for g0, g1 in self.grid) + len(self.prior_residuals(self.x))
        S = lil_matrix((nr, self.nparams), dtype=np.int8)
        # reprojection rows
        for o in range(self.nobs):
            k, jj = self.o_cam[o], self.o_pt[o]
            rows = [2 * o, 2 * o + 1]
            sl = self.cam_slices.get(k, {})
            cols = []
            for key in ("w", "dC", "s", "d"):
                if key in sl:
                    cols.extend(sl[key].tolist())
            if self.rho_idx is not None:
                cols.append(self.rho_idx)
            g0, g1 = self.grid[jj]
            n = g1 - g0 + 1
            if sl.get("d") is not None or self.rho_idx is not None:
                lo = self.o_tb[o] - self.d_hi[k] - g0
                hi = self.o_tb[o] - self.d_lo[k] - g0
            else:
                lo = hi = self.o_tb[o] - self.d0[k] - g0
            i_lo = max(1, int(np.floor(lo)) - 1)
            i_hi = min(n - 3, int(np.floor(hi)) + 1)
            for i in range(i_lo - 1, i_hi + 3):
                if 0 <= i < n:
                    base = self.traj_off[jj] + 3 * i
                    cols.extend([base, base + 1, base + 2])
            for r_ in rows:
                for c_ in cols:
                    S[r_, c_] = 1
        r0 = 2 * self.nobs
        for jj, (g0, g1) in enumerate(self.grid):
            n = g1 - g0 + 1
            for i in range(1, n - 1):
                for a in range(3):
                    for di in (-1, 0, 1):
                        S[r0 + a, self.traj_off[jj] + 3 * (i + di) + a] = 1
                r0 += 3
        # priors
        for k, sl in self.cam_slices.items():
            for key in ("w", "dC", "s", "d"):
                if key in sl:
                    for c_ in sl[key]:
                        S[r0, c_] = 1; r0 += 1
        if self.rho_idx is not None:
            S[r0, self.rho_idx] = 1; r0 += 1
        if self.gbasis is not None:
            for _ in range(self.gbasis.shape[0]):
                for c_ in self.gcols:
                    S[r0, c_] = 1
                r0 += 1
        assert r0 == nr, (r0, nr)
        return S.tocsr()

    def bounds(self):
        lo = np.full(self.nparams, -np.inf)
        hi = np.full(self.nparams, np.inf)
        for k, sl in self.cam_slices.items():
            if "d" in sl:
                lo[sl["d"]] = min(self.d0[k], self.d_prior[k]) - self.sync_bound
                hi[sl["d"]] = max(self.d0[k], self.d_prior[k]) + self.sync_bound
            if "w" in sl:
                lo[sl["w"]] = -0.2; hi[sl["w"]] = 0.2
                lo[sl["dC"]] = -1.0; hi[sl["dC"]] = 1.0
            if "s" in sl:
                lo[sl["s"]] = -0.1; hi[sl["s"]] = 0.1
        if self.rho_idx is not None:
            lo[self.rho_idx] = -1.5; hi[self.rho_idx] = 1.5
        return lo, hi

    def _irls_weights(self, x):
        r = self.reproj_residuals(x)
        a = np.linalg.norm(r, axis=1)
        self.w_obs = np.where(a <= self.hub_r, 1.0, self.hub_r / np.maximum(a, 1e-12))
        ra = self.acc_residuals(x)
        aa = np.linalg.norm(ra, axis=1)
        self.w_acc = np.where(aa <= self.hub_a, 1.0, self.hub_a / np.maximum(aa, 1e-12))

    def solve(self, irls_iters=4, max_nfev=60, verbose=0):
        S = self.sparsity()
        lo, hi = self.bounds()
        x = np.clip(self.x, lo + 1e-9, hi - 1e-9)
        self.w_obs = np.ones(self.nobs); self.w_acc = None
        hist = []
        for it in range(irls_iters):
            self._irls_weights(x)
            res = least_squares(self.fun, x, jac_sparsity=S, bounds=(lo, hi), method="trf",
                                x_scale="jac", max_nfev=max_nfev, verbose=verbose, tr_solver="lsmr")
            x = res.x
            hist.append(float(res.cost))
        self.x = x
        self.history = hist
        return self

    def robust_cost(self, x=None):
        """Huber objective (independent of the current IRLS weights)."""
        x = self.x if x is None else x
        def hub(a, c):
            return np.where(a <= c, 0.5 * a ** 2, c * a - 0.5 * c ** 2)
        a = np.linalg.norm(self.reproj_residuals(x), axis=1)
        aa = np.linalg.norm(self.acc_residuals(x), axis=1)
        pr = self.prior_residuals(x)
        return float(hub(a, self.hub_r).sum() + hub(aa, self.hub_a).sum() + 0.5 * np.sum(pr ** 2))

    # ------------------------------------------------------------------
    def report(self):
        r = self.reproj_residuals() * self.o_sig[:, None]
        e = np.linalg.norm(r, axis=1)
        cams = self.cameras()
        out = dict(ref_cam=self.ds.cam_ids[self.ref], rho=float(self.rho()),
                   d=dict(zip(self.ds.cam_ids, np.round(self.offsets_corr(), 3).tolist())), cams={})
        for k in range(self.C):
            m = self.o_cam == k
            row = dict(n=int(m.sum()), med_px=float(np.median(e[m])) if m.any() else np.nan)
            sl = self.cam_slices.get(k, {})
            if "w" in sl:
                row["rot_deg"] = float(np.rad2deg(np.linalg.norm(self.x[sl["w"]])))
                row["dC_cm"] = (100 * self.x[sl["dC"]]).round(1).tolist()
            if "s" in sl:
                row["focal_pct"] = float(100 * (np.exp(self.x[sl["s"]][0]) - 1))
            out["cams"][self.ds.cam_ids[k]] = row
        out["med_px_all"] = float(np.median(e))
        return out
