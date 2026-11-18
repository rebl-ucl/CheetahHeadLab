"""Joint (integer + sub-frame) camera synchronisation from manual labels.

For every camera pair (a, b) the symmetric point-to-epipolar-line distance of
corresponding labels is evaluated as a function of the relative time shift
s = d_b - d_a (camera b's track is linearly interpolated in time for fractional
s).  Pairwise curves do not depend on any third camera, so one badly synced view
cannot contaminate the estimate of another (the failure mode of a greedy
leave-one-out scan).  Global offsets are then the joint minimiser of the summed
pairwise costs: exhaustive over integers, then continuous coordinate descent.

Motion along an epipolar line is invisible to that pair; such pairs produce flat
curves and therefore contribute nothing — no explicit weighting needed.
"""
from __future__ import annotations

import itertools

import numpy as np

from .diagnostics import HEAD, TORSO


def fundamental(ca, cb):
    """F with x_b^T F x_a = 0 for undistorted MATLAB (1-based) pixels."""
    Rr = cb.R @ ca.R.T
    tr = cb.t - Rr @ ca.t
    tx = np.array([[0, -tr[2], tr[1]], [tr[2], 0, -tr[0]], [-tr[1], tr[0], 0]])
    Ka = ca.K.copy(); Ka[0, 2] += 1; Ka[1, 2] += 1
    Kb = cb.K.copy(); Kb[0, 2] += 1; Kb[1, 2] += 1
    return np.linalg.inv(Kb).T @ tx @ Rr @ np.linalg.inv(Ka)


def _interp_track(ud_b, rows):
    """Linear interpolation of a (N,P,2) track at fractional rows (M,) -> (M,P,2)."""
    N = ud_b.shape[0]
    i0 = np.floor(rows).astype(int)
    w = (rows - i0)[:, None, None]
    ok = (i0 >= 0) & (i0 + 1 < N)
    out = np.full((len(rows),) + ud_b.shape[1:], np.nan)
    i0c = np.clip(i0, 0, N - 2)
    a = ud_b[i0c]
    b = ud_b[i0c + 1]
    v = (1 - w) * a + w * b
    exact = np.isclose(w[:, 0, 0], 0.0)
    v[exact] = a[exact]  # integer shift: no need for the next frame
    out[ok] = v[ok]
    return out


class PairCost:
    """Epipolar cost of one camera pair as a function of relative shift."""

    def __init__(self, ds, a, b, uds, sigs, joints, cap=10.0, cams=None, row_filter=None):
        cams = ds.cams if cams is None else cams
        self.a, self.b = a, b
        self.F = fundamental(cams[a], cams[b])
        self.base = ds.offsets[b] - ds.offsets[a]
        self.ua = uds[a][:, joints].copy()
        self.ub = uds[b][:, joints].copy()
        if row_filter is not None:
            tb = np.arange(self.ua.shape[0]) + 1 - ds.offsets[a]
            self.ua[~row_filter(a, tb)] = np.nan
            tbb = np.arange(self.ub.shape[0]) + 1 - ds.offsets[b]
            self.ub[~row_filter(b, tbb)] = np.nan
        self.sa = sigs[a][:, joints]
        self.sb = sigs[b][:, joints]
        self.cap = cap
        self.rows_a = np.arange(self.ua.shape[0])

    def per_obs(self, s):
        """Per-observation capped cost (Na,P) at relative shift s (NaN = no pair)."""
        rows_b = self.rows_a + self.base + s
        ub = _interp_track(self.ub, rows_b.astype(float))
        sb = _interp_track(self.sb[..., None].repeat(2, -1), rows_b.astype(float))[..., 0]
        ua = self.ua
        m = np.all(np.isfinite(ua), -1) & np.all(np.isfinite(ub), -1)
        d = np.full(m.shape, np.nan)
        if m.any():
            xa = np.concatenate([ua[m], np.ones((m.sum(), 1))], 1)
            xb = np.concatenate([ub[m], np.ones((m.sum(), 1))], 1)
            lb = xa @ self.F.T
            la = xb @ self.F
            num = np.abs(np.sum(xb * lb, 1))
            dist = 0.5 * (num / np.hypot(lb[:, 0], lb[:, 1]) + num / np.hypot(la[:, 0], la[:, 1]))
            sig = np.sqrt(self.sa[m] ** 2 + sb[m] ** 2)
            d[m] = np.minimum(dist / sig, self.cap)
        return d

    def curve(self, shifts, common=False, min_n=15):
        P = [self.per_obs(s) for s in shifts]
        if common:
            m = np.logical_and.reduce([np.isfinite(p) for p in P])
            vals = np.array([p[m].mean() if m.sum() >= min_n else np.nan for p in P])
            n = np.full(len(shifts), int(m.sum()))
        else:
            vals = np.array([np.nanmean(p) if np.isfinite(p).sum() >= min_n else np.nan for p in P])
            n = np.array([int(np.isfinite(p).sum()) for p in P])
        return vals, n


def estimate_offsets(ds, R=5, joints=None, cams=None, fine_step=0.05, ref=None, verbose=False, row_filter=None):
    """Return dict with integer and fractional offset corrections (frames).

    corrections are ADDED to ds.offsets (camera row = abs frame + offset).
    """
    cams = ds.cams if cams is None else cams
    if joints is None:
        joints = [ds.pidx(p) for p in HEAD + TORSO + ["l_shoulder", "r_shoulder"] if ds.pidx(p) is not None]
    _, uds, sigs = ds.camera_tracks(cams)
    C = ds.ncams
    def _count(k):
        v = np.isfinite(uds[k][:, joints, 0])
        if row_filter is not None:
            tb = np.arange(v.shape[0]) + 1 - ds.offsets[k]
            v = v & row_filter(k, tb)[:, None]
        return int(v.sum())
    has = [_count(k) >= 30 for k in range(C)]
    active = [k for k in range(C) if has[k]]
    if ref is None:
        ref = max(active, key=_count)
    pairs = {}
    for a, b in itertools.combinations(active, 2):
        pairs[(a, b)] = PairCost(ds, a, b, uds, sigs, joints, cams=cams, row_filter=row_filter)

    # ---- stage A: integer curves over +-2R (available obs), exhaustive search ----
    ish = np.arange(-2 * R, 2 * R + 1)
    curvesA = {}
    for key, pc in pairs.items():
        v, n = pc.curve(ish, common=False)
        w = np.nanmedian(n)
        curvesA[key] = (v, n)
    others = [k for k in active if k != ref]
    grid = np.arange(-R, R + 1)
    combos = np.array(list(itertools.product(grid, repeat=len(others)))) if others else np.zeros((1, 0), int)
    off = np.zeros((len(combos), C), int)
    for i, k in enumerate(others):
        off[:, k] = combos[:, i]
    total = np.zeros(len(combos))
    for (a, b), (v, n) in curvesA.items():
        d = off[:, b] - off[:, a]
        vv = v[d + 2 * R]
        nn = n[d + 2 * R]
        # missing pair values: neutral (pair median)
        fill = np.nanmedian(v) if np.isfinite(v).any() else 0.0
        vv = np.where(np.isfinite(vv), vv, fill)
        total += vv * np.nanmedian(n)
    best = off[np.argmin(total)].astype(float)

    # ---- stage B: continuous refinement on common-set curves near the solution ----
    fs = np.round(np.arange(-2.0, 2.0 + 1e-9, fine_step), 6)
    curvesB = {}
    for (a, b), pc in pairs.items():
        d0 = best[b] - best[a]
        v, n = pc.curve(d0 + fs, common=True)
        if not np.isfinite(v).any():
            v, n = pc.curve(d0 + fs, common=False)
        curvesB[(a, b)] = (d0, v, n)

    def tot_cost(o):
        t = 0.0
        for (a, b), (d0, v, n) in curvesB.items():
            x = (o[b] - o[a]) - d0
            if abs(x) > 2.0 or not np.isfinite(v).any():
                t += (np.nanmax(v) if np.isfinite(v).any() else 0.0) * n[0]
                continue
            vi = np.interp(x, fs, np.where(np.isfinite(v), v, np.nanmax(v)))
            t += vi * n[0]
        return t

    o = best.copy()
    cur = tot_cost(o)
    for _ in range(30):
        changed = False
        for k in others:
            cand = o[k] + np.arange(-1.0, 1.0 + 1e-9, fine_step)
            vals = []
            for c in cand:
                oo = o.copy(); oo[k] = c
                vals.append(tot_cost(oo))
            i = int(np.argmin(vals))
            if vals[i] < cur - 1e-9 and abs(cand[i] - o[k]) > 1e-9:
                o[k] = cand[i]; cur = vals[i]; changed = True
        if not changed:
            break
    # per-camera curvature (confidence): cost increase for +-0.5 frame
    conf = {}
    for k in others:
        oo1 = o.copy(); oo1[k] += 0.5
        oo2 = o.copy(); oo2[k] -= 0.5
        conf[ds.cam_ids[k]] = float(0.5 * (tot_cost(oo1) + tot_cost(oo2)) - cur)
    out = dict(ref=ds.cam_ids[ref], active=[ds.cam_ids[k] for k in active],
               int_corr={ds.cam_ids[k]: int(best[k]) for k in range(C)},
               frac_corr={ds.cam_ids[k]: float(o[k]) for k in range(C)},
               cost_int=float(total.min()), cost_zero=float(total[np.all(off == 0, 1)][0]),
               confidence=conf,
               pair_curves={f"{ds.cam_ids[a]}-{ds.cam_ids[b]}": dict(shifts=ish.tolist(), cost=np.round(v, 3).tolist(),
                                                                   n=n.tolist())
                            for (a, b), (v, n) in curvesA.items()})
    if verbose:
        print("ref cam", out["ref"], "int", out["int_corr"], "frac", {k: round(v, 2) for k, v in out["frac_corr"].items()})
        print("cost zero %.1f -> int %.1f" % (out["cost_zero"], out["cost_int"]), "conf", {k: round(v, 1) for k, v in conf.items()})
    return out


def aligned_arrays(ds, corr, cams=None):
    """Undistorted labels on the common frame grid ds.frames after adding the
    (possibly fractional) per-camera corrections `corr` (frames) to ds.offsets.
    Fractional shifts use linear interpolation of the camera's own track.
    Returns ud (C,F,P,2) and sigma (C,F,P)."""
    cams = ds.cams if cams is None else cams
    _, uds, sigs = ds.camera_tracks(cams)
    C, F, P = ds.ncams, len(ds.frames), len(ds.points)
    ud = np.full((C, F, P, 2), np.nan)
    sg = np.full((C, F, P), np.nan)
    for k in range(C):
        rows = ds.frames + ds.offsets[k] + corr[k] - 1.0
        ud[k] = _interp_track(uds[k], rows)
        s = _interp_track(np.repeat(sigs[k][..., None], 2, -1), rows)[..., 0]
        sg[k] = s
    return ud, sg
