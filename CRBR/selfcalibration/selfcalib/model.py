"""Python replica of CRBR/reconstruction/cheetah3D.m and head/torso axis angles from pose vectors."""
from __future__ import annotations

import numpy as np

BODY = ["l_eye", "r_eye", "nose", "l_ear", "r_ear", "spine", "neck_base"]


def Rz(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1.0]])


def Rx(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[1.0, 0, 0], [0, c, s], [0, -s, c]])


def Ry(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])


def cheetah3d(p, sH, sT):
    """p: 12-vector pose (MATLAB order). Returns (R_head, R_torso, dict of points)."""
    tx, ty, tz, x, y, z, tx2, ty2, tz2, x2, y2, z2 = p
    RH = Rz(tz) @ Rx(tx) @ Ry(ty)
    RT = Rz(tz2) @ Rx(tx2) @ Ry(ty2)
    ph = np.array([x, y, z]); pt = np.array([x2, y2, z2])
    pts = {
        "l_eye": ph + RH @ [0, sH * 0.03, 0],
        "r_eye": ph + RH @ [0, -sH * 0.03, 0],
        "l_ear": ph + RH @ [sH * -0.0725, sH * 0.0608, sH * 0.0441],
        "r_ear": ph + RH @ [sH * -0.0725, sH * -0.0608, sH * 0.0441],
        "nose": ph + RH @ [sH * 0.055, 0, -sH * 0.055],
        "neck_base": pt + RT @ [sT * 0.04, 0, sT * 0.10],
        "spine": pt + RT @ [-sT * 0.33, 0, sT * 0.10],
    }
    return RH, RT, pts


def head_torso_angles(xbar):
    """World yaw/pitch (deg) of head and torso forward (+X local) axes. xbar (F,12)."""
    out = np.full((len(xbar), 4), np.nan)
    for i, p in enumerate(xbar):
        if not np.all(np.isfinite(p)):
            continue
        RH, RT, _ = cheetah3d(p, 1.0, 1.0)
        for k, R in enumerate((RH, RT)):
            f = R @ np.array([1.0, 0, 0])
            out[i, 2 * k] = np.degrees(np.arctan2(f[1], f[0]))
            out[i, 2 * k + 1] = np.degrees(np.arctan2(f[2], np.hypot(f[0], f[1])))
    return out  # columns: head yaw, head pitch, torso yaw, torso pitch
