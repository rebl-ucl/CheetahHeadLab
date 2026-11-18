"""Write a self-calibrated dataset for the subsequent CRBR stages.

Output folder contents
  AcinoSet_data/<n>_cam_scene_sba.json   refined K/D/R/t (same JSON layout)
  AcinoSet_data/offset.xlsx              integer offsets  off + round(d)
  cam<N>_labels_wide.csv                 labels time-aligned for the sub-frame
                                         remainder phi = d - round(d) (and the
                                         rolling-shutter term, if fitted)
  cam<N>_llh_wide.csv, limits.xlsx       copied
  cam<N>_mjpeg.avi                       copied (overlay figures and videos)
  selfcalib_report.json                  every estimated parameter

Time alignment: row r of camera c is captured at t = (r+1-off'_c) - phi_c
[+ rho*(v/H-0.5)] but the pipeline reads it as frame f = r+1-off'_c.  Each label
is moved by the image displacement the fitted 3-D motion predicts between t and
f:  u' = u + pi_c(X(f)) - pi_c(X(t)).  Only the predicted motion is removed; the
label's own residual (annotation noise) is preserved.  |phi| <= 0.5 frame, so
the correction never exceeds half the inter-frame image motion.
"""
from __future__ import annotations

import json
import os
import shutil

import numpy as np
import pandas as pd

from .camera import save_scene_json
from .dataset import BODY


def _traj_shift(b, jj, cam, t_from, t_to):
    g0, g1 = b.grid[jj]
    ok = (t_from >= g0 + 1) & (t_from <= g1 - 2) & (t_to >= g0 + 1) & (t_to <= g1 - 2)
    sh = np.full((len(t_from), 2), np.nan)
    if ok.any():
        Xa = b.eval_points(b.x, t_from[ok], np.full(ok.sum(), jj))
        Xb = b.eval_points(b.x, t_to[ok], np.full(ok.sum(), jj))
        sh[ok] = cam.project(Xb, distort=True) - cam.project(Xa, distort=True)
    return sh


def export_dataset(ds, fits, out_folder, report=None, copy_videos=True):
    """fits: list of fitted STBA objects sharing cameras/d/rho (body fit first,
    optional lure fit second). Returns dict summary of applied shifts."""
    body = fits[0]
    cams = body.cameras()
    d = body.offsets_corr().copy()
    rho = body.rho()
    # anchor the time gauge on the reference camera (keeps limits.xlsx meaning)
    d -= np.round(d[body.ref])
    d_int = np.round(d).astype(int)
    phi = d - d_int
    new_off = ds.offsets + d_int

    os.makedirs(os.path.join(out_folder, "AcinoSet_data"), exist_ok=True)
    save_scene_json(os.path.join(out_folder, "AcinoSet_data", f"{ds.ncams}_cam_scene_sba.json"), cams,
                    template_js={k: v for k, v in ds.scene_js.items() if k != "cameras"},
                    extra={"selfcalib": "refined by CRBR spatiotemporal self-calibration; original: " + ds.json_path})
    off_tbl = pd.DataFrame({"cam": ds.cam_ids, "offset": new_off})
    off_tbl.to_excel(os.path.join(out_folder, "AcinoSet_data", "offset.xlsx"), index=False)
    for extra in ("FTE_start_frame.xlsx",):
        src = os.path.join(ds.folder, "AcinoSet_data", extra)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(out_folder, "AcinoSet_data", extra))
    if os.path.isfile(os.path.join(ds.folder, "limits.xlsx")):
        shutil.copy2(os.path.join(ds.folder, "limits.xlsx"), os.path.join(out_folder, "limits.xlsx"))

    stats = {}
    for k, cid in enumerate(ds.cam_ids):
        lt = ds.label_tables[k].copy()
        ht = ds.llh_tables[k].copy()
        N = len(lt)
        r = np.arange(N)
        f_new = r + 1 - new_off[k]
        shifts_all = []
        for p in list(ds.points):
            if f"{p}_x" not in lt.columns:
                continue
            u = np.stack([pd.to_numeric(lt[f"{p}_x"], errors="coerce").values,
                          pd.to_numeric(lt[f"{p}_y"], errors="coerce").values], 1).astype(float)
            fin = np.all(np.isfinite(u), 1)
            t_cap = f_new - phi[k] + rho * (u[:, 1] / ds.H - 0.5)
            sh = np.full_like(u, np.nan)
            fit = next((b for b in fits if p in b.points), None)
            if fit is not None and fin.any():
                jj = fit.points.index(p)
                idx = np.where(fin)[0]
                sh[idx] = _traj_shift(fit, jj, cams[k], t_cap[idx], f_new[idx].astype(float))
            # fallback: camera's own central-difference image velocity
            need = fin & ~np.all(np.isfinite(sh), 1)
            if need.any():
                vel = np.full_like(u, np.nan)
                vel[1:-1] = (u[2:] - u[:-2]) / 2.0
                fb = need & np.all(np.isfinite(vel), 1)
                sh[fb] = vel[fb] * (f_new[fb] - t_cap[fb])[:, None]
                still = need & ~fb
                sh[still] = 0.0
            sh[~fin] = 0.0
            unew = u + sh
            lt[f"{p}_x"] = np.where(fin, unew[:, 0], lt[f"{p}_x"])
            lt[f"{p}_y"] = np.where(fin, unew[:, 1], lt[f"{p}_y"])
            shifts_all.append(np.linalg.norm(sh[fin], axis=1))
        # pad missing body columns so the MATLAB loader's asserts pass
        for p in BODY:
            if f"{p}_x" not in lt.columns:
                lt[f"{p}_x"] = np.nan
                lt[f"{p}_y"] = np.nan
            if f"{p}_LLH" not in ht.columns:
                ht[f"{p}_LLH"] = 0.0
        lt.to_csv(os.path.join(out_folder, f"cam{cid}_labels_wide.csv"), index=False)
        ht.to_csv(os.path.join(out_folder, f"cam{cid}_llh_wide.csv"), index=False)
        v = np.concatenate(shifts_all) if shifts_all else np.zeros(0)
        stats[cid] = dict(phi=float(phi[k]), int_offset=int(new_off[k]), orig_offset=int(ds.offsets[k]),
                          median_shift_px=float(np.median(v)) if v.size else 0.0,
                          p95_shift_px=float(np.percentile(v, 95)) if v.size else 0.0)
        if copy_videos:
            for ext in ("_mjpeg.avi",):
                src = os.path.join(ds.folder, f"cam{cid}{ext}")
                dst = os.path.join(out_folder, f"cam{cid}{ext}")
                if os.path.isfile(src) and not os.path.exists(dst):
                    shutil.copy2(src, dst)
    rep = dict(source=ds.folder, cam_ids=ds.cam_ids, orig_offsets=ds.offsets.tolist(),
               new_offsets=new_off.tolist(), d=d.tolist(), phi=phi.tolist(), rho=float(rho),
               fit=body.report(), label_shift=stats)
    if report:
        rep.update(report)
    with open(os.path.join(out_folder, "selfcalib_report.json"), "w") as f:
        json.dump(rep, f, indent=1, default=float)
    return rep


def export_baseline(ds, out_folder, copy_videos=True):
    """Copy of the original dataset with only the missing-column padding, so the
    current pipeline can be run on it as a like-for-like 'before'."""
    os.makedirs(os.path.join(out_folder, "AcinoSet_data"), exist_ok=True)
    for f in os.listdir(os.path.join(ds.folder, "AcinoSet_data")):
        if f.endswith(".json") or f.endswith(".xlsx"):
            shutil.copy2(os.path.join(ds.folder, "AcinoSet_data", f), os.path.join(out_folder, "AcinoSet_data", f))
    if os.path.isfile(os.path.join(ds.folder, "limits.xlsx")):
        shutil.copy2(os.path.join(ds.folder, "limits.xlsx"), os.path.join(out_folder, "limits.xlsx"))
    for k, cid in enumerate(ds.cam_ids):
        lt = ds.label_tables[k].copy()
        ht = ds.llh_tables[k].copy()
        for p in BODY:
            if f"{p}_x" not in lt.columns:
                lt[f"{p}_x"] = np.nan
                lt[f"{p}_y"] = np.nan
            if f"{p}_LLH" not in ht.columns:
                ht[f"{p}_LLH"] = 0.0
        lt.to_csv(os.path.join(out_folder, f"cam{cid}_labels_wide.csv"), index=False)
        ht.to_csv(os.path.join(out_folder, f"cam{cid}_llh_wide.csv"), index=False)
        if copy_videos:
            src = os.path.join(ds.folder, f"cam{cid}_mjpeg.avi")
            dst = os.path.join(out_folder, f"cam{cid}_mjpeg.avi")
            if os.path.isfile(src) and not os.path.exists(dst):
                shutil.copy2(src, dst)
