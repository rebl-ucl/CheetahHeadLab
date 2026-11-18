"""Fast per-frame robust rigid head/torso fit (Python replica of the CRBR
framewise stage) used for camera-deletion stability tests.

Mirrors cheetahError.m: undistorted pinhole residuals / sigma (x3 when
LLH < 0.85), per-landmark cross-camera consensus weights (>=3 cams, scale floor
2, cut-off 4x scale, redescending) and a Huber weight (c = 3).  Head and torso
are independent in the model, so they are fitted separately; torso roll about
the spine-neck axis is unobservable and irrelevant to the forward axis.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .camera import rodrigues
from .triangulate import triangulate

TEMPLATE_HEAD = {"l_eye": [0, 0.03, 0], "r_eye": [0, -0.03, 0], "nose": [0.055, 0, -0.055],
                 "l_ear": [-0.0725, 0.0608, 0.0441], "r_ear": [-0.0725, -0.0608, 0.0441]}
TEMPLATE_TORSO = {"neck_base": [0.04, 0, 0.10], "spine": [-0.33, 0, 0.10]}


def estimate_scales(ds, X):
    def d(a, b):
        ia, ib = ds.pidx(a), ds.pidx(b)
        if ia is None or ib is None:
            return np.array([])
        v = np.linalg.norm(X[:, ia] - X[:, ib], axis=-1)
        return v[np.isfinite(v)]
    de = d("l_eye", "r_eye")
    il, ir, inn = ds.pidx("l_eye"), ds.pidx("r_eye"), ds.pidx("nose")
    dn = np.linalg.norm(X[:, inn] - 0.5 * (X[:, il] + X[:, ir]), axis=-1)
    dn = dn[np.isfinite(dn)]
    s_eye = np.median(de) / 0.06 if de.size else np.nan
    s_nose = np.median(dn) / np.hypot(0.055, 0.055) if dn.size else np.nan
    sH = np.nanmean([s_eye, s_nose])
    if np.isfinite(s_eye) and np.isfinite(s_nose) and abs(s_eye - s_nose) / s_nose > 0.2:
        sH = s_nose
    dt = d("spine", "neck_base")
    sT = np.median(dt) / 0.37 if dt.size else 1.0
    return float(sH), float(sT)


def _kabsch(A, B):
    """R, t minimising |R A + t - B| (rows are points)."""
    ca, cb = A.mean(0), B.mean(0)
    H = (A - ca).T @ (B - cb)
    U, _, Vt = np.linalg.svd(H)
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, cb - R @ ca


class FrameFitter:
    def __init__(self, ds, cams, ud, sigma, scale_head, scale_torso, cons_C=4.0, cons_floor=2.0, huber=3.0):
        self.ds, self.cams, self.ud, self.sigma = ds, cams, ud, sigma
        self.sH, self.sT = scale_head, scale_torso
        self.cC, self.cF, self.hub = cons_C, cons_floor, huber
        self.P = [np.hstack([c.R, c.t[:, None]]) for c in cams]

    def _residual(self, x, R0, names, Q, f, cams_ok):
        R = rodrigues(x[:3]) @ R0
        o = x[3:6]
        raw = []
        for jn, (n, q) in enumerate(zip(names, Q)):
            Xw = o + R @ q
            j = self.ds.pidx(n)
            for k in cams_ok[jn]:
                c = self.cams[k]
                Xc = c.R @ Xw + c.t
                p = np.array([c.K[0, 0] * Xc[0] / Xc[2] + c.K[0, 2] + 1, c.K[1, 1] * Xc[1] / Xc[2] + c.K[1, 2] + 1])
                raw.append(((p - self.ud[k, f, j]) / self.sigma[k, f, j], jn, k))
        if not raw:
            return np.zeros(1)
        mags = np.array([np.linalg.norm(r[0]) for r in raw])
        w = np.ones(len(raw))
        jns = np.array([r[1] for r in raw])
        for jn in np.unique(jns):
            m = jns == jn
            if m.sum() >= 3:
                sc = max(np.median(mags[m]), self.cF)
                ratio = mags[m] / sc
                w[m] = np.where(ratio <= self.cC, (1 - (ratio / self.cC) ** 2) ** 2, 0.0)
        hw = np.where(mags <= self.hub, 1.0, self.hub / np.maximum(mags, 1e-12))
        out = np.concatenate([r[0] * np.sqrt(wi * hi) for r, wi, hi in zip(raw, w, hw)])
        return out

    def fit_segment(self, f, segment, exclude_cams=(), R_init=None, o_init=None):
        tmpl = TEMPLATE_HEAD if segment == "head" else TEMPLATE_TORSO
        s = self.sH if segment == "head" else self.sT
        names = [n for n in tmpl if self.ds.pidx(n) is not None]
        Q = [s * np.array(tmpl[n], float) for n in names]
        cams_ok = []
        Xtri = []
        for n in names:
            j = self.ds.pidx(n)
            ks = [k for k in range(self.ds.ncams) if k not in exclude_cams and np.all(np.isfinite(self.ud[k, f, j]))]
            cams_ok.append(ks)
            if len(ks) >= 2:
                X, _, dep = triangulate([self.cams[k] for k in ks], self.ud[ks, f, j], w=1.0 / self.sigma[ks, f, j])
                Xtri.append(X if np.all(dep > 0) else np.full(3, np.nan))
            else:
                Xtri.append(np.full(3, np.nan))
        nobs = sum(len(k) for k in cams_ok)
        Xtri = np.array(Xtri)
        okp = np.all(np.isfinite(Xtri), 1)
        if segment == "head":
            if okp.sum() < 3 or nobs < 6:
                return None
            R0, o0 = _kabsch(np.array(Q)[okp], Xtri[okp])
        else:
            if okp.sum() < 2 or nobs < 4:
                return None
            v = Xtri[names.index("neck_base")] - Xtri[names.index("spine")]
            v = v / np.linalg.norm(v)
            # level, zero-roll frame whose +X maps onto v
            yaw = np.arctan2(v[1], v[0]); pitch = np.arctan2(v[2], np.hypot(v[0], v[1]))
            cz, sz, cp, sp = np.cos(yaw), np.sin(yaw), np.cos(pitch), np.sin(pitch)
            R0 = np.array([[cz * cp, -sz, -cz * sp], [sz * cp, cz, -sz * sp], [sp, 0, cp]])
            mid = 0.5 * (Xtri[0] + Xtri[1])
            o0 = mid - R0 @ (0.5 * (Q[0] + Q[1]))
        if R_init is not None:
            R0, o0 = R_init, o_init
        x0 = np.concatenate([np.zeros(3), o0])
        res = least_squares(self._residual, x0, args=(R0, names, Q, f, cams_ok), method="lm" if nobs * 2 >= 6 else "trf",
                            xtol=1e-10, ftol=1e-10, max_nfev=200)
        R = rodrigues(res.x[:3]) @ R0
        fwd = R @ np.array([1.0, 0, 0])
        return dict(R=R, o=res.x[3:6], fwd=fwd, nobs=nobs, cost=float(res.cost))


def yaw_pitch(fwd):
    return np.degrees(np.arctan2(fwd[1], fwd[0])), np.degrees(np.arctan2(fwd[2], np.hypot(fwd[0], fwd[1])))


def angle_between(a, b):
    return np.degrees(np.arccos(np.clip(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1)))


def deletion_test(ds, cams, ud, sigma, frames_mask, scale_head, scale_torso):
    """Per camera deletion: angle (deg) between forward axes with/without the camera."""
    ff = FrameFitter(ds, cams, ud, sigma, scale_head, scale_torso)
    out = {"head": {}, "torso": {}}
    base = {}
    for seg in ("head", "torso"):
        for f in np.where(frames_mask)[0]:
            r = ff.fit_segment(f, seg)
            if r is not None:
                base[(seg, f)] = r
    for k, cid in enumerate(ds.cam_ids):
        for seg in ("head", "torso"):
            angs, dyaw, dpitch = [], [], []
            for f in np.where(frames_mask)[0]:
                if (seg, f) not in base:
                    continue
                b = base[(seg, f)]
                names = list((TEMPLATE_HEAD if seg == "head" else TEMPLATE_TORSO).keys())
                if not any(ds.pidx(n) is not None and np.all(np.isfinite(ud[k, f, ds.pidx(n)])) for n in names):
                    continue  # camera contributes nothing this frame
                r = ff.fit_segment(f, seg, exclude_cams=(k,), R_init=b["R"], o_init=b["o"])
                if r is None:
                    continue
                angs.append(angle_between(b["fwd"], r["fwd"]))
                y0, p0 = yaw_pitch(b["fwd"]); y1, p1 = yaw_pitch(r["fwd"])
                dyaw.append(abs((y1 - y0 + 180) % 360 - 180)); dpitch.append(abs(p1 - p0))
            out[seg][cid] = dict(n=len(angs), med_angle=float(np.median(angs)) if angs else np.nan,
                                 p95_angle=float(np.percentile(angs, 95)) if angs else np.nan,
                                 med_dyaw=float(np.median(dyaw)) if dyaw else np.nan,
                                 med_dpitch=float(np.median(dpitch)) if dpitch else np.nan)
    return out, base
