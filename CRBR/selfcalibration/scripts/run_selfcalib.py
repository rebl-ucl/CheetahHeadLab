"""Self-calibrate one dataset (sync + camera extrinsics) and export a corrected
copy for the subsequent CRBR reconstruction stages.

usage:
  run_selfcalib.py SRC OUT [BASELINE_OUT] [--cv] [--sync-only] [--tight]

  SRC           original dataset folder (never modified)
  OUT           corrected copy (refined JSON, corrected offsets, time-aligned labels)
  BASELINE_OUT  optional like-for-like copy of the original (only missing ear
                columns padded) for before/after comparison with the pipeline
  --cv          run temporal-block cross-validation (original / sync-only /
                sync+calibration) and keep calibration refinement only if it
                generalises better than sync-only on held-out frames+cameras
  --sync-only   never refine the calibration
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np  # noqa: E402

from selfcalib.dataset import Dataset  # noqa: E402
from selfcalib.diagnostics import swap_check  # noqa: E402
from selfcalib.evaluate import cv_evaluate, summarize  # noqa: E402
from selfcalib.export import export_baseline, export_dataset  # noqa: E402
from selfcalib.pipeline import fit_geometry, method_for_cv  # noqa: E402
from selfcalib.stba import STBA  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("out"); ap.add_argument("baseline", nargs="?")
ap.add_argument("--fps", type=float, default=None, help="recording frame rate from manifest.json")
ap.add_argument("--cv", action="store_true")
ap.add_argument("--sync-only", action="store_true")
ap.add_argument("--tight", action="store_true", help="sigma_rot 0.5 deg, sigma_pos 5 cm")
ap.add_argument("--min-obs", type=int, default=0)
ap.add_argument("--window", type=int, default=None,
                help="fit only labels within this many frames of limits.xlsx (non-stationary sync)")
a = ap.parse_args()
if os.path.exists(a.out):
    ap.error("Output already exists. Choose a new folder to protect original inputs and prior results.")
if a.baseline and os.path.exists(a.baseline):
    ap.error("Baseline output already exists. Choose a new folder.")

t0 = time.time()
ds = Dataset(a.src)
fps = a.fps if a.fps is not None else (90.0 if ds.W == 1920 else 120.0)
print(ds.summary(), flush=True)
cfg = dict(calib_min_obs=a.min_obs)
if a.tight:
    cfg.update(sigma_rot_deg=0.5, sigma_pos_m=0.05)
report = dict(fps=fps, cfg=cfg, window=a.window)
win_filter = None
if a.window is not None and ds.limits is not None:
    lo_w, hi_w = ds.limits[0] - a.window, ds.limits[1] + a.window
    win_filter = lambda k, p, tb: (tb >= lo_w) & (tb <= hi_w)
    print(f"fitting only frames {lo_w}..{hi_w}", flush=True)

# ---- optional held-out model check ------------------------------------------------
use_calib = not a.sync_only
if a.cv:
    cvres = {}
    for mode in ("orig", "sync", "sync_calib"):
        base_m = method_for_cv(cfg, mode)
        m_ = base_m if win_filter is None else (
            lambda d_, f_, base_m=base_m: base_m(d_, lambda k, p, tb: f_(k, p, tb) & win_filter(k, p, tb)))
        r, folds = cv_evaluate(ds, m_, K=4, fps=fps,
                               t_range=None if win_filter is None else (lo_w, hi_w))
        cvres[mode] = dict(all=summarize(r), head=summarize(r, ["l_eye", "r_eye", "nose", "l_ear", "r_ear"]),
                           folds=folds)
        s = cvres[mode]["all"]
        print(f"CV {mode:10s}: held-out median {s['all'][0]:.2f}px mean {s['all_mean'][0]:.2f}px (n={s['all'][1]}) "
              + " ".join(f"c{c}:{v[0]:.1f}" for c, v in s.items() if c not in ("all", "all_mean")), flush=True)
    sc, so, sg = cvres["sync_calib"]["all"], cvres["sync"]["all"], cvres["orig"]["all"]
    if not a.sync_only:
        use_calib = (sc["all"][0] <= so["all"][0]) and (sc["all_mean"][0] <= 1.1 * so["all_mean"][0])
    # per-camera safeguard: a camera whose refined calibration predicts its own
    # held-out labels clearly worse than sync-only keeps its original calibration
    freeze = []
    if use_calib:
        for k, cid in enumerate(ds.cam_ids):
            a_, b_ = sc.get(cid, (np.nan, 0)), so.get(cid, (np.nan, 0))
            if b_[1] >= 30 and np.isfinite(a_[0]) and np.isfinite(b_[0]) and a_[0] > 1.25 * b_[0]:
                freeze.append(k)
    if freeze:
        cfg["freeze_calib_cams"] = freeze
        report["frozen_calibration_cams"] = [ds.cam_ids[k] for k in freeze]
        print("calibration kept original for cams", [ds.cam_ids[k] for k in freeze], flush=True)
        base_m = method_for_cv(cfg, "sync_calib")
        m_ = base_m if win_filter is None else (
            lambda d_, f_, base_m=base_m: base_m(d_, lambda k, p, tb: f_(k, p, tb) & win_filter(k, p, tb)))
        r, folds = cv_evaluate(ds, m_, K=4, fps=fps, t_range=None if win_filter is None else (lo_w, hi_w))
        cvres["sync_calib_frozen"] = dict(all=summarize(r), head=summarize(r, ["l_eye", "r_eye", "nose", "l_ear", "r_ear"]),
                                          folds=folds)
        s_ = cvres["sync_calib_frozen"]["all"]
        print(f"CV sync_calib (frozen {report['frozen_calibration_cams']}): held-out median {s_['all'][0]:.2f}px "
              f"mean {s_['all_mean'][0]:.2f}px " + " ".join(f"c{c}:{v[0]:.1f}" for c, v in s_.items() if c not in ("all", "all_mean")), flush=True)
        # keep the freeze only if it also wins overall: camera calibrations are
        # coupled, so freezing one view can push the error into another
        if s_["all"][0] <= sc["all"][0] and s_["all_mean"][0] <= sc["all_mean"][0]:
            sc = s_
        else:
            cfg.pop("freeze_calib_cams")
            report["frozen_calibration_cams_rejected"] = report.pop("frozen_calibration_cams")
            print("freeze rejected: it does not improve overall held-out error", flush=True)
    chosen = sc if use_calib else so
    keep_original = (sg["all"][0] < chosen["all"][0]) and (sg["all_mean"][0] < chosen["all_mean"][0])
    report["cv"] = cvres
    report["cv_choice"] = "original" if keep_original else ("sync_calib" if use_calib else "sync_only")
    if keep_original:
        print("WARNING: original geometry predicts held-out labels best; exporting it unchanged", flush=True)
    print("CV choice:", report["cv_choice"], flush=True)

# ---- final fit on all data -----------------------------------------------------------
s2, s1, log = fit_geometry(ds, obs_filter=win_filter, cfg=dict(cfg, refine_calib=use_calib), fps=fps, verbose=True)
body = s2 if s2 is not None else s1
if report.get("cv_choice") == "original":
    # identity geometry: original cameras and offsets (structure-only fit for the export)
    body = STBA(ds, refine_calib=False, refine_sync=False, fps=fps, obs_filter=win_filter, sync_bound=0.0)
    body.solve(irls_iters=2, max_nfev=30)
report["fit_log"] = log
print("final geometry:", json.dumps(body.report(), default=float), flush=True)

# ---- lure trajectory (geometry fixed; sync only for cameras without body labels) ----
fits = [body]
if "lure" in ds.points:
    no_body = [k for k in range(ds.ncams) if not body.cam_has[k]]
    lure = STBA(ds, points=["lure"], cams=body.cameras(), d0=body.offsets_corr(), refine_calib=False,
                refine_sync=bool(no_body), sync_cams=no_body, fps=fps, sync_bound=4.0,
                acc_ms2=300.0, huber_reproj=2.0, ref_cam=body.ref, obs_filter=win_filter)
    lure.solve(irls_iters=3, max_nfev=40)
    if no_body:
        d = body.offsets_corr()
        for k in no_body:
            d[k] = lure.offsets_corr()[k]
        body.d0 = d
        print("lure-based sync for cams", [ds.cam_ids[k] for k in no_body], np.round(d, 2).tolist())
    report["lure_fit"] = lure.report()
    fits.append(lure)

rep = export_dataset(ds, fits, a.out, report=report, copy_videos=False)
print("exported", a.out, "offsets", ds.offsets.tolist(), "->", rep["new_offsets"],
      "phi", np.round(rep["phi"], 2).tolist())
for c, s in rep["label_shift"].items():
    print(f"  cam{c}: offset {s['orig_offset']}->{s['int_offset']} phi={s['phi']:+.2f} "
          f"label shift median {s['median_shift_px']:.1f}px p95 {s['p95_shift_px']:.1f}px")

# ---- post-correction left/right swap screen (flag only; labels are not altered) ----
dsc = Dataset(a.out)
sw = swap_check(dsc)
flag = [s for s in sw if s["e_swap"] < 0.5 * s["e_keep"] and s["e_keep"] - s["e_swap"] > 15]
rep["swap_screen"] = dict(tested=len(sw), flagged=flag)
with open(os.path.join(a.out, "selfcalib_report.json"), "w") as f:
    json.dump(rep, f, indent=1, default=float)
print(f"swap screen: {len(flag)} flagged of {len(sw)} tested" +
      ("".join(f"\n  {x['pair']} cam{x['cam']} frame {x['frame']}: {x['e_keep']:.0f}px -> {x['e_swap']:.0f}px if swapped"
               for x in flag[:20])))
if a.baseline:
    export_baseline(ds, a.baseline, copy_videos=False)
    print("baseline exported", a.baseline)
print(f"done in {time.time() - t0:.0f}s")
