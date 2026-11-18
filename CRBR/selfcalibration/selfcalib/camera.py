"""Camera model matching the MATLAB reconstruction pipeline exactly.

AcinoSet scene JSONs store OpenCV fisheye (Kannala-Brandt) intrinsics K, D and a
world->camera extrinsic (R, t):  X_c = R X_w + t.

MATLAB's ``cameraIntrinsicsFromOpenCV`` turns a 4-element D into a
``cameraIntrinsicsKB`` object whose principal point is the OpenCV one + 1
(MATLAB pixels are 1-based).  Verified numerically against R2026a:
``undistortPoints`` == ``cv2.fisheye.undistortPoints(uv-1, K, D, P=K) + 1`` and
``world2img`` == pinhole projection with (cx+1, cy+1).  Every function here
works in that MATLAB 1-based pixel convention, so labels can be used as stored.
"""
from __future__ import annotations

import numpy as np


def rodrigues(w):
    """Rotation matrix from rotation vector(s). w: (3,) -> (3,3)."""
    w = np.asarray(w, float)
    th = np.linalg.norm(w)
    if th < 1e-12:
        W = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
        return np.eye(3) + W
    k = w / th
    Kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * Kx + (1 - np.cos(th)) * (Kx @ Kx)


def log_so3(R):
    """Rotation vector of a rotation matrix."""
    c = np.clip((np.trace(R) - 1) / 2, -1.0, 1.0)
    th = np.arccos(c)
    if th < 1e-9:
        return np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) / 2
    if np.pi - th < 1e-6:
        # near pi: use the symmetric part
        A = (R + np.eye(3)) / 2
        k = np.sqrt(np.clip(np.diag(A), 0, None))
        i = int(np.argmax(k))
        k[i] = np.sqrt(A[i, i])
        for j in range(3):
            if j != i:
                k[j] = A[i, j] / k[i]
        return th * k / np.linalg.norm(k)
    return th / (2 * np.sin(th)) * np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])


def kb_distort(xn, D):
    """Kannala-Brandt distortion of normalised pinhole coords. xn: (...,2)."""
    D = np.asarray(D, float).reshape(4)
    r = np.sqrt(np.sum(xn ** 2, axis=-1))
    th = np.arctan(r)
    th2 = th * th
    thd = th * (1 + D[0] * th2 + D[1] * th2 ** 2 + D[2] * th2 ** 3 + D[3] * th2 ** 4)
    with np.errstate(invalid="ignore", divide="ignore"):
        s = np.where(r > 1e-12, thd / np.where(r > 1e-12, r, 1.0), 1.0)
    return xn * s[..., None]


def kb_undistort(xd, D, iters=30):
    """Inverse of kb_distort (Newton on theta), matching cv2.fisheye.undistortPoints."""
    D = np.asarray(D, float).reshape(4)
    thd = np.sqrt(np.sum(xd ** 2, axis=-1))
    th = thd.copy()
    for _ in range(iters):
        th2 = th * th
        f = th * (1 + D[0] * th2 + D[1] * th2 ** 2 + D[2] * th2 ** 3 + D[3] * th2 ** 4) - thd
        fp = 1 + 3 * D[0] * th2 + 5 * D[1] * th2 ** 2 + 7 * D[2] * th2 ** 3 + 9 * D[3] * th2 ** 4
        th = th - f / fp
    with np.errstate(invalid="ignore", divide="ignore"):
        s = np.where(thd > 1e-12, np.tan(th) / np.where(thd > 1e-12, thd, 1.0), 1.0)
    return xd * s[..., None]


class Camera:
    """One calibrated camera (OpenCV K/D, world->camera R/t), MATLAB pixel convention."""

    def __init__(self, K, D, R, t, size):
        self.K = np.asarray(K, float).copy()
        self.D = np.asarray(D, float).reshape(4).copy()
        self.R = np.asarray(R, float).copy()
        self.t = np.asarray(t, float).reshape(3).copy()
        self.size = tuple(int(s) for s in size)  # (W, H)

    def copy(self):
        return Camera(self.K, self.D, self.R, self.t, self.size)

    @property
    def center(self):
        return -self.R.T @ self.t

    # ---- pixel <-> normalised (MATLAB 1-based) ----
    def pix_to_norm(self, uv):
        K = self.K
        x = (uv[..., 0] - 1 - K[0, 2] - K[0, 1] * 0) / K[0, 0]
        y = (uv[..., 1] - 1 - K[1, 2]) / K[1, 1]
        return np.stack([x, y], axis=-1)

    def norm_to_pix(self, xn):
        K = self.K
        u = K[0, 0] * xn[..., 0] + K[0, 2] + 1
        v = K[1, 1] * xn[..., 1] + K[1, 2] + 1
        return np.stack([u, v], axis=-1)

    def undistort(self, uv_raw):
        """Raw (distorted) MATLAB pixels -> undistorted MATLAB pixels (same K)."""
        return self.norm_to_pix(kb_undistort(self.pix_to_norm(uv_raw), self.D))

    def to_cam(self, X):
        return X @ self.R.T + self.t

    def project(self, X, distort=True):
        """World points (...,3) -> MATLAB pixels; raw/distorted if distort else undistorted."""
        Xc = self.to_cam(X)
        xn = Xc[..., :2] / Xc[..., 2:3]
        if distort:
            xn = kb_distort(xn, self.D)
        return self.norm_to_pix(xn)

    def depth(self, X):
        return self.to_cam(X)[..., 2]

    def ray_norm(self, uv_undist):
        """Normalised pinhole coordinates of undistorted pixels."""
        return self.pix_to_norm(uv_undist)


def load_scene_json(path):
    import json
    with open(path) as f:
        js = json.load(f)
    W, H = js["camera_resolution"]
    cams = []
    for c in js["cameras"]:
        cams.append(Camera(np.array(c["k"], float), np.array(c["d"], float).reshape(-1),
                           np.array(c["r"], float), np.array(c["t"], float).reshape(-1), (W, H)))
    return cams, js


def save_scene_json(path, cams, template_js=None, extra=None):
    import json
    js = dict(template_js) if template_js is not None else {}
    W, H = cams[0].size
    js["camera_resolution"] = [W, H]
    out = []
    for c in cams:
        out.append({
            "k": c.K.tolist(),
            "d": [[float(v)] for v in c.D],
            "r": c.R.tolist(),
            "t": [[float(v)] for v in c.t],
        })
    js["cameras"] = out
    if extra:
        js.update(extra)
    with open(path, "w") as f:
        json.dump(js, f, indent=2)
