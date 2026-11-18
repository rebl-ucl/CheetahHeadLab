"""Multiview triangulation in undistorted MATLAB pixels (DLT + Gauss-Newton)."""
from __future__ import annotations

import numpy as np


def _pinhole_proj_jac(cam, X):
    """Undistorted pixel projection of one point and its 2x3 Jacobian w.r.t. X."""
    Xc = cam.R @ X + cam.t
    z = Xc[2]
    fx, fy = cam.K[0, 0], cam.K[1, 1]
    u = fx * Xc[0] / z + cam.K[0, 2] + 1
    v = fy * Xc[1] / z + cam.K[1, 2] + 1
    J = np.array([[fx / z, 0, -fx * Xc[0] / z ** 2],
                  [0, fy / z, -fy * Xc[1] / z ** 2]]) @ cam.R
    return np.array([u, v]), J, z


def triangulate(cams, uv, w=None, refine=True, iters=10):
    """Triangulate one point.

    cams : list of Camera (length M >= 2)
    uv   : (M,2) undistorted MATLAB pixels
    w    : (M,) weights (1/sigma), optional
    Returns X (3,), per-view residual vectors (M,2) in px, and per-view depths.
    """
    M = len(cams)
    if w is None:
        w = np.ones(M)
    A = np.zeros((2 * M, 4))
    for i, cam in enumerate(cams):
        xn = cam.pix_to_norm(uv[i])
        P = np.hstack([cam.R, cam.t[:, None]])
        A[2 * i] = w[i] * (xn[0] * P[2] - P[0])
        A[2 * i + 1] = w[i] * (xn[1] * P[2] - P[1])
    _, _, Vt = np.linalg.svd(A)
    Xh = Vt[-1]
    if abs(Xh[3]) < 1e-12:
        return np.full(3, np.nan), np.full((M, 2), np.nan), np.full(M, np.nan)
    X = Xh[:3] / Xh[3]
    if refine:
        for _ in range(iters):
            r = []
            J = []
            for i, cam in enumerate(cams):
                p, Ji, _ = _pinhole_proj_jac(cam, X)
                r.append(w[i] * (p - uv[i]))
                J.append(w[i] * Ji)
            r = np.concatenate(r)
            J = np.vstack(J)
            try:
                dx = -np.linalg.lstsq(J, r, rcond=None)[0]
            except np.linalg.LinAlgError:
                break
            X = X + dx
            if np.linalg.norm(dx) < 1e-7:
                break
    res = np.zeros((M, 2))
    dep = np.zeros(M)
    for i, cam in enumerate(cams):
        p, _, z = _pinhole_proj_jac(cam, X)
        res[i] = p - uv[i]
        dep[i] = z
    return X, res, dep


def triangulate_dataset(ds, cams=None, ud=None, cam_mask=None, min_views=2):
    """Triangulate every (frame, point) of a Dataset from all visible cameras.

    Returns X (F,P,3), res (C,F,P,2) in-sample residuals, nviews (F,P).
    """
    cams = ds.cams if cams is None else cams
    ud = ds.ud if ud is None else ud
    C, F, P, _ = ud.shape
    X = np.full((F, P, 3), np.nan)
    res = np.full((C, F, P, 2), np.nan)
    nv = np.zeros((F, P), int)
    for f in range(F):
        for j in range(P):
            ks = [k for k in range(C) if np.all(np.isfinite(ud[k, f, j])) and (cam_mask is None or cam_mask[k, f])]
            nv[f, j] = len(ks)
            if len(ks) < min_views:
                continue
            Xj, r, dep = triangulate([cams[k] for k in ks], ud[ks, f, j], w=1.0 / ds.sigma[ks, f, j])
            if not np.all(np.isfinite(Xj)) or np.any(dep <= 0):
                continue
            X[f, j] = Xj
            res[ks, f, j] = r
    return X, res, nv


def loo_residuals(ds, cams=None, ud=None, min_views=3):
    """Leave-one-camera-out residuals: triangulate from the OTHER views, reproject.

    Returns loo (C,F,P,2) residual vectors (projection - label), undistorted px,
    and the held-out camera depth of the LOO point zloo (C,F,P).
    """
    cams = ds.cams if cams is None else cams
    ud = ds.ud if ud is None else ud
    C, F, P, _ = ud.shape
    loo = np.full((C, F, P, 2), np.nan)
    zloo = np.full((C, F, P), np.nan)
    for f in range(F):
        for j in range(P):
            ks = [k for k in range(C) if np.all(np.isfinite(ud[k, f, j]))]
            if len(ks) < min_views:
                continue
            for k in ks:
                others = [m for m in ks if m != k]
                Xj, _, dep = triangulate([cams[m] for m in others], ud[others, f, j],
                                         w=1.0 / ds.sigma[others, f, j])
                if not np.all(np.isfinite(Xj)) or np.any(dep <= 0):
                    continue
                p = cams[k].project(Xj, distort=False)
                loo[k, f, j] = p - ud[k, f, j]
                zloo[k, f, j] = cams[k].depth(Xj)
    return loo, zloo
