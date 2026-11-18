"""The self-calibration recipe shared by the production tool and the CV harness.

    stage 1  sync only: multi-start (original offsets / pairwise-epipolar
             estimate), joint multiview objective, pick the lowest robust cost
    stage 2  joint sync + per-camera extrinsic refinement, similarity gauge
             fixed, cameras with too little data keep their calibration
"""
from __future__ import annotations

import numpy as np

from .stba import STBA
from .sync import estimate_offsets

DEFAULTS = dict(
    stage1_bound=4.0,        # frames reachable from each start in stage 1
    stage2_bound=1.0,        # frames reachable around the stage-1 estimate
    sigma_sync=2.0,          # sync prior (frames), centred on the original offsets
    sigma_rot_deg=0.5,       # extrinsic rotation prior (CV-selected, see README)
    sigma_pos_m=0.05,        # camera-centre prior (CV-selected)
    calib_min_obs=0,         # cameras with fewer observations keep their calibration
    freeze_calib_cams=(),    # camera indices that keep their original calibration
    gauge="similarity",      # remove common-mode (7-DOF similarity) corrections
    refine_calib=True,
    acc_ms2=100.0,           # smoothness prior on landmark trajectories
)


def fit_geometry(ds, obs_filter=None, cfg=None, fps=None, verbose=False):
    """Returns (stage2 STBA or None, stage1 STBA, log dict)."""
    c = dict(DEFAULTS)
    if cfg:
        c.update(cfg)
    fps = fps or (90.0 if ds.W == 1920 else 120.0)
    log = dict(cfg=c)
    starts = {"original": np.zeros(ds.ncams)}
    try:
        rf = None if obs_filter is None else (lambda k, tb: obs_filter(k, None, tb))
        r = estimate_offsets(ds, row_filter=rf)
        starts["epipolar"] = np.array([r["frac_corr"][cid] for cid in ds.cam_ids])
    except Exception as e:  # pragma: no cover - logged, not fatal
        log["epipolar_error"] = repr(e)
    best = None
    log["stage1"] = {}
    for name, d0 in starts.items():
        b = STBA(ds, d0=d0, refine_calib=False, refine_sync=True, fps=fps, obs_filter=obs_filter,
                 sync_bound=c["stage1_bound"], sigma_sync=c["sigma_sync"], acc_ms2=c["acc_ms2"])
        b.solve(irls_iters=3, max_nfev=40)
        cost = b.robust_cost()
        log["stage1"][name] = dict(cost=cost, d=b.offsets_corr().tolist())
        if verbose:
            print(f"  stage1 start={name}: cost {cost:.1f} d={np.round(b.offsets_corr(), 2).tolist()}", flush=True)
        if best is None or cost < best[0]:
            best = (cost, b, name)
    s1 = best[1]
    log["stage1_choice"] = best[2]
    if not c["refine_calib"]:
        return None, s1, log
    s2 = STBA(ds, d0=s1.offsets_corr(), refine_calib=True, refine_sync=True, fps=fps, obs_filter=obs_filter,
              sync_bound=c["stage2_bound"], sigma_sync=c["sigma_sync"], sigma_rot_deg=c["sigma_rot_deg"],
              sigma_pos_m=c["sigma_pos_m"], gauge=c["gauge"], calib_min_obs=c["calib_min_obs"],
              freeze_calib_cams=c["freeze_calib_cams"],
              acc_ms2=c["acc_ms2"])
    s2.solve(irls_iters=4, max_nfev=60)
    log["stage2"] = s2.report()
    return s2, s1, log


def method_for_cv(cfg=None, mode="sync_calib"):
    """Adapter for evaluate.cv_evaluate: returns (cams, d, rho)."""
    def m(ds, filt):
        if mode == "orig":
            return ds.cams, np.zeros(ds.ncams), 0.0
        s2, s1, _ = fit_geometry(ds, obs_filter=filt, cfg=dict(cfg or {}, refine_calib=(mode == "sync_calib")))
        if s2 is None:
            return ds.cams, s1.offsets_corr(), 0.0
        return s2.cameras(), s2.offsets_corr(), s2.rho()
    return m
