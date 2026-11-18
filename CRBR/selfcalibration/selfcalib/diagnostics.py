"""Geometric diagnostics that separate calibration, synchronisation and label errors.

All residuals are leave-one-camera-out (LOO): the landmark is triangulated from
the OTHER cameras and projected into the held-out camera, so a biased camera
cannot hide its own error by pulling the triangulation (paper note S2).

Signatures used to attribute the error:
  * sync:  r_LOO ~ -delta * v_img  (residual parallel to the camera's own image
           velocity, proportional to speed)  ->  fitted sub-frame offset delta
  * integer offset: LOO error minimised at a shift s != 0
  * calibration: residual NOT explained by velocity, structured in image
    position / depth (perpendicular component)
  * label swap: swapping a bilateral pair in one camera reduces its LOO error
"""
from __future__ import annotations

import numpy as np

from .triangulate import triangulate, loo_residuals

HEAD = ["l_eye", "r_eye", "nose", "l_ear", "r_ear"]
TORSO = ["spine", "neck_base"]


def _shift(arr, s):
    """arr (F,...) -> out[f] = arr[f+s] (NaN padded)."""
    out = np.full_like(arr, np.nan)
    F = arr.shape[0]
    if s >= 0:
        out[: F - s] = arr[s:]
    else:
        out[-s:] = arr[: F + s]
    return out


def loo_for_camera(ds, k, ud_k, cams=None, joints=None, min_other=2):
    """LOO residual of camera k given its (possibly shifted) labels ud_k (F,P,2)."""
    cams = ds.cams if cams is None else cams
    C, F, P, _ = ds.ud.shape
    joints = range(P) if joints is None else joints
    out = np.full((F, P, 2), np.nan)
    for f in range(F):
        for j in joints:
            if not np.all(np.isfinite(ud_k[f, j])):
                continue
            others = [m for m in range(C) if m != k and np.all(np.isfinite(ds.ud[m, f, j]))]
            if len(others) < min_other:
                continue
            X, _, dep = triangulate([cams[m] for m in others], ds.ud[others, f, j], w=1.0 / ds.sigma[others, f, j])
            if not np.all(np.isfinite(X)) or np.any(dep <= 0) or cams[k].depth(X) <= 0:
                continue
            out[f, j] = cams[k].project(X, distort=False) - ud_k[f, j]
    return out


def integer_offset_scan(ds, shifts=range(-3, 4), joints=None):
    """Per camera: median LOO error vs integer shift of that camera's labels."""
    C = ds.ncams
    joints = [ds.pidx(p) for p in HEAD + TORSO if ds.pidx(p) is not None] if joints is None else joints
    table = {}
    for k in range(C):
        per_s = {}
        masks = {}
        for s in shifts:
            r = loo_for_camera(ds, k, _shift(ds.ud[k], s), joints=joints)
            per_s[s] = r
            masks[s] = np.all(np.isfinite(r), axis=-1)
        common = np.logical_and.reduce([masks[s] for s in shifts])
        row = {}
        for s in shifts:
            m = common if common.sum() >= 10 else masks[s]
            e = np.linalg.norm(per_s[s][m], axis=-1)
            row[int(s)] = (float(np.median(e)) if e.size else np.nan, int(e.size))
        table[ds.cam_ids[k]] = row
    return table


def image_velocity(ds):
    """Central-difference image velocity (undistorted px/frame) per camera/frame/point."""
    v = np.full_like(ds.ud, np.nan)
    v[:, 1:-1] = (ds.ud[:, 2:] - ds.ud[:, :-2]) / 2.0
    return v


def sync_signature(ds, loo, vel, joints=None, min_speed=2.0):
    """Fit r_LOO = -delta * v per camera. Returns dict cam -> stats."""
    out = {}
    joints = list(range(ds.ud.shape[2])) if joints is None else joints
    for k, cid in enumerate(ds.cam_ids):
        r = loo[k][:, joints].reshape(-1, 2)
        v = vel[k][:, joints].reshape(-1, 2)
        m = np.all(np.isfinite(r), 1) & np.all(np.isfinite(v), 1)
        sp = np.linalg.norm(v, axis=1)
        m &= sp > min_speed
        if m.sum() < 8:
            out[cid] = dict(n=int(m.sum()))
            continue
        r, v, sp = r[m], v[m], sp[m]
        d_ls = -np.sum(r * v) / np.sum(v * v)
        # robust: IRLS Huber on scalar model
        d = d_ls
        for _ in range(20):
            e = r + d * v
            en = np.linalg.norm(e, axis=1)
            c = 1.345 * np.median(en) / 0.6745 + 1e-9
            w = np.where(en <= c, 1.0, c / en)
            d = -np.sum(w[:, None] * r * v) / np.sum(w[:, None] * v * v)
        par = np.sum(r * v, 1) / sp  # residual component along velocity (px)
        perp = (r[:, 0] * v[:, 1] - r[:, 1] * v[:, 0]) / sp
        tot = np.sum(r ** 2)
        expl = 1 - np.sum((r + d * v) ** 2) / tot if tot > 0 else np.nan
        out[cid] = dict(n=int(m.sum()), delta_ls=float(d_ls), delta_robust=float(d),
                        frac_var_explained=float(expl),
                        med_abs_par=float(np.median(np.abs(par))), med_abs_perp=float(np.median(np.abs(perp))),
                        med_speed=float(np.median(sp)))
    return out


def swap_check(ds, pairs=(("l_eye", "r_eye"), ("l_ear", "r_ear"), ("l_shoulder", "r_shoulder"))):
    """For each camera/frame where both members of a bilateral pair are labelled,
    compare LOO error of the labels as given vs swapped (other cams unchanged)."""
    res = []
    C, F, P, _ = ds.ud.shape
    for a, b in pairs:
        ia, ib = ds.pidx(a), ds.pidx(b)
        if ia is None or ib is None:
            continue
        for k in range(C):
            for f in range(F):
                ua, ub = ds.ud[k, f, ia], ds.ud[k, f, ib]
                if not (np.all(np.isfinite(ua)) and np.all(np.isfinite(ub))):
                    continue
                Xs = {}
                for j in (ia, ib):
                    others = [m for m in range(C) if m != k and np.all(np.isfinite(ds.ud[m, f, j]))]
                    if len(others) < 2:
                        break
                    X, _, dep = triangulate([ds.cams[m] for m in others], ds.ud[others, f, j])
                    if not np.all(np.isfinite(X)) or np.any(dep <= 0):
                        break
                    Xs[j] = ds.cams[k].project(X, distort=False)
                if len(Xs) < 2:
                    continue
                e_keep = np.linalg.norm(Xs[ia] - ua) + np.linalg.norm(Xs[ib] - ub)
                e_swap = np.linalg.norm(Xs[ia] - ub) + np.linalg.norm(Xs[ib] - ua)
                res.append(dict(pair=f"{a}/{b}", cam=ds.cam_ids[k], frame=int(ds.frames[f]),
                                e_keep=float(e_keep), e_swap=float(e_swap)))
    return res


def anatomical_distances(ds, X):
    """Distances between triangulated landmarks (m), per frame."""
    def d(a, b):
        ia, ib = ds.pidx(a), ds.pidx(b)
        if ia is None or ib is None:
            return None
        return np.linalg.norm(X[:, ia] - X[:, ib], axis=-1)
    out = {}
    for a, b in [("l_eye", "r_eye"), ("l_ear", "r_ear"), ("l_eye", "nose"), ("r_eye", "nose"),
                 ("l_eye", "l_ear"), ("r_eye", "r_ear"), ("spine", "neck_base"), ("nose", "neck_base")]:
        v = d(a, b)
        if v is None:
            continue
        v = v[np.isfinite(v)]
        if v.size < 3:
            continue
        q = np.percentile(v, [10, 25, 50, 75, 90])
        out[f"{a}-{b}"] = dict(n=int(v.size), p10=float(q[0]), p25=float(q[1]), median=float(q[2]),
                              p75=float(q[3]), p90=float(q[4]), rcv=float((q[3] - q[1]) / 1.349 / q[2]))
    return out


def per_camera_loo_stats(ds, loo, zloo, frame_mask=None):
    """Median LOO error (px) per camera & landmark group; depth relation."""
    out = {}
    groups = {"head": [ds.pidx(p) for p in HEAD if ds.pidx(p) is not None],
              "torso": [ds.pidx(p) for p in TORSO if ds.pidx(p) is not None],
              "shoulder": [ds.pidx(p) for p in ("l_shoulder", "r_shoulder") if ds.pidx(p) is not None],
              "lure": [ds.pidx("lure")] if ds.pidx("lure") is not None else []}
    for k, cid in enumerate(ds.cam_ids):
        row = {}
        for g, js in groups.items():
            if not js:
                continue
            r = loo[k][:, js] if frame_mask is None else loo[k][frame_mask][:, js]
            e = np.linalg.norm(r.reshape(-1, 2), axis=1)
            e = e[np.isfinite(e)]
            row[g] = (float(np.median(e)) if e.size else np.nan, int(e.size))
        js = groups["head"] + groups["torso"]
        z = zloo[k][:, js].ravel()
        z = z[np.isfinite(z)]
        row["median_depth_m"] = float(np.median(z)) if z.size else np.nan
        row["min_depth_m"] = float(np.min(z)) if z.size else np.nan
        out[cid] = row
    return out


def epipolar_pairs(ds, joints=None, frame_mask=None):
    """Symmetric point-to-epipolar-line distance (undistorted px) per camera pair."""
    C = ds.ncams
    joints = [ds.pidx(p) for p in HEAD + TORSO if ds.pidx(p) is not None] if joints is None else joints
    out = {}
    for a in range(C):
        for b in range(a + 1, C):
            ca, cb = ds.cams[a], ds.cams[b]
            Rr = cb.R @ ca.R.T
            tr = cb.t - Rr @ ca.t
            tx = np.array([[0, -tr[2], tr[1]], [tr[2], 0, -tr[0]], [-tr[1], tr[0], 0]])
            E = tx @ Rr
            Ka = ca.K.copy(); Ka[0, 2] += 1; Ka[1, 2] += 1
            Kb = cb.K.copy(); Kb[0, 2] += 1; Kb[1, 2] += 1
            Fm = np.linalg.inv(Kb).T @ E @ np.linalg.inv(Ka)
            ua = ds.ud[a][:, joints].reshape(-1, 2)
            ub = ds.ud[b][:, joints].reshape(-1, 2)
            if frame_mask is not None:
                ua = ds.ud[a][frame_mask][:, joints].reshape(-1, 2)
                ub = ds.ud[b][frame_mask][:, joints].reshape(-1, 2)
            m = np.all(np.isfinite(ua), 1) & np.all(np.isfinite(ub), 1)
            if m.sum() < 5:
                continue
            ha = np.c_[ua[m], np.ones(m.sum())]
            hb = np.c_[ub[m], np.ones(m.sum())]
            lb = ha @ Fm.T
            la = hb @ Fm
            db = np.abs(np.sum(hb * lb, 1)) / np.hypot(lb[:, 0], lb[:, 1])
            da = np.abs(np.sum(ha * la, 1)) / np.hypot(la[:, 0], la[:, 1])
            out[f"{ds.cam_ids[a]}-{ds.cam_ids[b]}"] = dict(n=int(m.sum()), med_b=float(np.median(db)),
                                                          med_a=float(np.median(da)))
    return out


def diagnose(ds, do_offsets=True, do_swaps=True):
    from .triangulate import triangulate_dataset
    X, res, nv = triangulate_dataset(ds)
    loo, zloo = loo_residuals(ds)
    vel = image_velocity(ds)
    lm = ds.limit_mask()
    body_j = [ds.pidx(p) for p in HEAD + TORSO if ds.pidx(p) is not None]
    rep = dict(name=ds.name, cams=ds.cam_ids, offsets=ds.offsets.tolist(), limits=ds.limits,
               frames=[int(ds.frames[0]), int(ds.frames[-1])], res=[ds.W, ds.H],
               missing_points=ds.missing_points)
    rep["loo_all"] = per_camera_loo_stats(ds, loo, zloo)
    rep["loo_limits"] = per_camera_loo_stats(ds, loo, zloo, frame_mask=lm)
    rep["sync_signature_body"] = sync_signature(ds, loo, vel, joints=body_j)
    rep["epipolar_limits"] = epipolar_pairs(ds, frame_mask=lm)
    rep["anatomy"] = anatomical_distances(ds, X)
    # views per frame within limits
    rep["views_limits"] = {p: float(np.mean(nv[lm, j])) for j, p in enumerate(ds.points)}
    if do_offsets:
        rep["offset_scan"] = integer_offset_scan(ds)
    if do_swaps:
        sw = swap_check(ds)
        flagged = [s for s in sw if s["e_swap"] < 0.5 * s["e_keep"] and s["e_keep"] - s["e_swap"] > 15]
        rep["swap_tested"] = len(sw)
        rep["swap_flagged"] = flagged
    return rep, dict(X=X, res=res, nv=nv, loo=loo, zloo=zloo, vel=vel)
