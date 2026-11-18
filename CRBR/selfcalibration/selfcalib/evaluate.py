"""Honest out-of-sample evaluation: temporal-block CV + held-out-camera prediction.

For each of K contiguous time blocks:
  1. fit the method (cameras / sync) on observations OUTSIDE the block
     (plus a guard gap so the smoothness prior cannot leak information);
  2. inside the block, for every camera c: fit the landmark trajectories from
     the OTHER cameras only (camera parameters fixed at the training estimate),
     then predict camera c's labels at its own capture times.
The prediction error (raw px) therefore uses neither the held-out frames nor
the held-out camera's labels in those frames.
"""
from __future__ import annotations

import numpy as np

from .stba import STBA


def obs_time_range(ds, points):
    raws, _, _ = ds.camera_tracks()
    ts = []
    for k in range(ds.ncams):
        for p in points:
            j = ds.pidx(p)
            if j is None:
                continue
            r = np.where(np.all(np.isfinite(raws[k][:, j]), -1))[0]
            ts.append(r + 1 - ds.offsets[k])
    t = np.concatenate(ts)
    return t.min(), t.max()


def predict_block(ds, cams, d, rho, block, test_cam, points, fps, min_support=2):
    """Structure from cameras != test_cam inside `block`; predict test_cam labels."""
    lo, hi = block

    def filt(k, p, tb):
        return (tb >= lo) & (tb <= hi)

    b = STBA(ds, points=points, cams=cams, d0=d, refine_calib=False, refine_sync=False,
             rolling_shutter=False, fps=fps, obs_filter=filt, drop_cams_obs=(test_cam,), sync_bound=0.0)
    if b.nobs < 10:
        return None
    b.solve(irls_iters=2, max_nfev=30)
    # test observations
    raws, _, sigs = ds.camera_tracks(cams)
    t_all = b.obs_times()
    errs = []
    for jj, p in enumerate(b.points):
        j = ds.pidx(p)
        r = np.where(np.all(np.isfinite(raws[test_cam][:, j]), -1))[0]
        tb = r + 1 - ds.offsets[test_cam]
        keep = (tb >= lo) & (tb <= hi)
        r, tb = r[keep], tb[keep]
        if len(r) == 0:
            continue
        u = raws[test_cam][r, j]
        t = tb - d[test_cam] + rho * (u[:, 1] / ds.H - 0.5)
        g0, g1 = b.grid[jj]
        if g1 < g0:
            continue
        inside = (t >= g0 + 2) & (t <= g1 - 2)
        # support: >= min_support other cameras observe the point within +-1 frame
        m = b.o_pt == jj
        sup = np.array([len(np.unique(b.o_cam[m][np.abs(t_all[m] - ti) <= 1.0])) for ti in t])
        ok = inside & (sup >= min_support)
        if not ok.any():
            continue
        X = b.eval_points(b.x, t[ok], np.full(ok.sum(), jj))
        pred = cams[test_cam].project(X, distort=True)
        e = np.linalg.norm(pred - u[ok], axis=1)
        for ee, ff in zip(e, tb[ok]):
            errs.append((jj, p, float(ee), int(ff)))
    return errs


def cv_evaluate(ds, method, K=4, guard=3, points=None, fps=90.0, verbose=False, scheme="contiguous",
                block_len=12, t_range=None):
    """method(ds, obs_filter) -> (cams, d, rho). Returns per-camera held-out errors.

    scheme='contiguous'  : K contiguous time blocks (extrapolation to unseen
                           positions; strict for cameras that see the animal briefly)
    scheme='interleaved' : blocks of block_len frames dealt cyclically to K folds
                           (interpolation; matches production, where calibration
                           is fitted on all frames including the analysis window)
    """
    from .stba import CALIB_POINTS
    points = [p for p in (CALIB_POINTS if points is None else points) if p in ds.points]
    tmin, tmax = obs_time_range(ds, points)
    if t_range is not None:
        tmin, tmax = max(tmin, t_range[0]), min(tmax, t_range[1])
    if scheme == "contiguous":
        e = np.linspace(tmin, tmax + 1, K + 1)
        fold_blocks = [[(e[k], e[k + 1] - 1e-9)] for k in range(K)]
    else:
        starts = np.arange(tmin, tmax + 1, block_len)
        fold_blocks = [[] for _ in range(K)]
        for i, b0 in enumerate(starts):
            fold_blocks[i % K].append((b0, b0 + block_len - 1e-9))
    results = {c: [] for c in ds.cam_ids}
    fold_params = []
    for kf in range(K):
        blocks = fold_blocks[kf]

        def train_filter(k, p, tb, blocks=blocks):
            keep = np.ones(len(tb), bool)
            for lo, hi in blocks:
                keep &= (tb < lo - guard) | (tb > hi + guard)
            return keep

        cams, d, rho = method(ds, train_filter)
        fold_params.append(dict(d=np.asarray(d).tolist(), rho=float(rho)))
        for lo, hi in blocks:
            for c in range(ds.ncams):
                errs = predict_block(ds, cams, np.asarray(d), rho, (lo, hi), c, points, fps)
                if errs:
                    results[ds.cam_ids[c]].extend(errs)
        if verbose:
            print(f"  fold {kf}: d={np.round(d, 2).tolist()} rho={rho:.2f}", flush=True)
    return results, fold_params


def summarize(results, points_filter=None):
    out = {}
    allv = []
    for c, errs in results.items():
        v = np.array([e[2] for e in errs if points_filter is None or e[1] in points_filter])
        out[c] = (float(np.median(v)) if v.size else np.nan, int(v.size))
        allv.append(v)
    allv = np.concatenate(allv) if allv else np.zeros(0)
    out["all"] = (float(np.median(allv)) if allv.size else np.nan, int(allv.size))
    out["all_mean"] = (float(np.mean(allv)) if allv.size else np.nan, int(allv.size))
    return out
