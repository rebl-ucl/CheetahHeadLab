"""Load a CheetahHeadLab dataset folder with the CRBR camera and annotation conventions.

Conventions mirrored from multiviewRecon.m / lureReconTriangulation.m:
  * camera IDs = sorted N of cam<N>_labels_wide.csv; calibration file
    AcinoSet_data/<ncams>_cam_scene_sba.json, camera k of the JSON <-> camIDs[k]
  * offsets: AcinoSet_data/offset.xlsx column 'offset', taken BY POSITION k
  * absolute frame f (1-based) of camera k lives in label row f + offset_k
    (1-based data row, i.e. python row f + offset_k - 1)
  * an observation is visible when LLH > 0.5, finite and both coords > 1;
    sigma is multiplied by 3 when LLH < 0.85
"""
from __future__ import annotations

import glob
import os
import re

import numpy as np
import pandas as pd

from .camera import load_scene_json

BODY = ["l_eye", "r_eye", "nose", "l_ear", "r_ear", "spine", "neck_base"]
EXTRA = ["l_shoulder", "r_shoulder"]
LURE = "lure"
ALL_POINTS = BODY + EXTRA + [LURE]

# a-priori label sigmas (px) — production values from the paper (eq. S27);
# shoulders/lure are not in the production model, given conservative values.
SIGMA = {"l_eye": 3.0, "r_eye": 3.0, "nose": 2.0, "l_ear": 5.0, "r_ear": 5.0,
         "spine": 6.0, "neck_base": 6.0, "l_shoulder": 8.0, "r_shoulder": 8.0, "lure": 6.0}

BILATERAL = [("l_eye", "r_eye"), ("l_ear", "r_ear"), ("l_shoulder", "r_shoulder")]


class Dataset:
    def __init__(self, folder, offsets_override=None, json_override=None):
        self.folder = os.path.abspath(folder)
        self.name = os.path.basename(self.folder.rstrip("/"))
        files = glob.glob(os.path.join(self.folder, "cam*_labels_wide.csv"))
        self.cam_ids = sorted(int(re.search(r"cam(\d+)_labels_wide", os.path.basename(f)).group(1)) for f in files)
        self.ncams = len(self.cam_ids)
        jpath = json_override or os.path.join(self.folder, "AcinoSet_data", f"{self.ncams}_cam_scene_sba.json")
        self.json_path = jpath
        self.cams, self.scene_js = load_scene_json(jpath)
        assert len(self.cams) == self.ncams, f"{jpath}: {len(self.cams)} cameras vs {self.ncams} label files"
        W, H = self.cams[0].size
        self.W, self.H = W, H

        # offsets (by position, as multiviewRecon does)
        opath = os.path.join(self.folder, "AcinoSet_data", "offset.xlsx")
        if offsets_override is not None:
            self.offsets = np.asarray(offsets_override, int)
        elif os.path.isfile(opath):
            ot = pd.read_excel(opath)
            self.offsets = np.asarray(ot["offset"].values[: self.ncams], float)
            assert np.all(np.isfinite(self.offsets)), f"non-finite offsets {self.offsets}"
            self.offsets = self.offsets.astype(int)
        else:
            self.offsets = np.zeros(self.ncams, int)

        # limits
        lpath = os.path.join(self.folder, "limits.xlsx")
        self.limits = None
        if os.path.isfile(lpath):
            lt = pd.read_excel(lpath)
            self.limits = (int(lt.iloc[0, 0]), int(lt.iloc[0, 1]))

        # labels
        self.label_tables = []
        self.llh_tables = []
        for c in self.cam_ids:
            lt = pd.read_csv(os.path.join(self.folder, f"cam{c}_labels_wide.csv"))
            ht = pd.read_csv(os.path.join(self.folder, f"cam{c}_llh_wide.csv"))
            self.label_tables.append(lt)
            self.llh_tables.append(ht)
        self.points = [p for p in ALL_POINTS if any(f"{p}_x" in t.columns for t in self.label_tables)]
        self.missing_points = [p for p in ALL_POINTS if p not in self.points]
        self._build_arrays()

    # ------------------------------------------------------------------
    def _build_arrays(self, offsets=None):
        off = self.offsets if offsets is None else np.asarray(offsets, int)
        nrows = np.array([len(t) for t in self.label_tables])
        f_min = int(max(1, np.max(1 - off)))
        f_max = int(np.min(nrows - off))
        self.frames = np.arange(f_min, f_max + 1)  # absolute frames
        F, C, P = len(self.frames), self.ncams, len(self.points)
        raw = np.full((C, F, P, 2), np.nan)
        llh = np.zeros((C, F, P))
        for k in range(C):
            lt, ht = self.label_tables[k], self.llh_tables[k]
            rows = self.frames + off[k] - 1
            for j, p in enumerate(self.points):
                if f"{p}_x" not in lt.columns:
                    continue
                xy = np.stack([pd.to_numeric(lt[f"{p}_x"], errors="coerce").values[rows],
                               pd.to_numeric(lt[f"{p}_y"], errors="coerce").values[rows]], axis=1).astype(float)
                lcol = f"{p}_LLH"
                l = pd.to_numeric(ht[lcol], errors="coerce").values[rows].astype(float) if lcol in ht.columns else np.zeros(F)
                raw[k, :, j] = xy
                llh[k, :, j] = np.nan_to_num(l, nan=0.0)
        vis = (llh > 0.5) & np.all(np.isfinite(raw), axis=-1) & np.all(raw > 1, axis=-1)
        raw[~vis] = np.nan
        self.raw = raw
        self.llh = llh
        self.vis = vis
        sig = np.array([SIGMA[p] for p in self.points])
        self.sigma = np.where(llh < 0.85, 3.0, 1.0) * sig[None, None, :]
        self.cur_offsets = off.copy()
        # undistorted (MATLAB convention, fixed calibration)
        ud = np.full_like(raw, np.nan)
        for k in range(C):
            m = vis[k]
            ud[k][m] = self.cams[k].undistort(raw[k][m])
        self.ud = ud

    def camera_tracks(self, cams=None):
        """Per-camera full-length label tracks indexed by the camera's own row.

        Returns lists (one per camera) of raw (N_k,P,2), undistorted (N_k,P,2) and
        sigma (N_k,P). Row r (0-based) of camera k shows absolute frame
        r + 1 - offset_k under the current offsets.
        """
        cams = self.cams if cams is None else cams
        raws, uds, sigs = [], [], []
        sig0 = np.array([SIGMA[p] for p in self.points])
        for k in range(self.ncams):
            lt, ht = self.label_tables[k], self.llh_tables[k]
            N = len(lt)
            raw = np.full((N, len(self.points), 2), np.nan)
            llh = np.zeros((N, len(self.points)))
            for j, p in enumerate(self.points):
                if f"{p}_x" not in lt.columns:
                    continue
                raw[:, j, 0] = pd.to_numeric(lt[f"{p}_x"], errors="coerce").values
                raw[:, j, 1] = pd.to_numeric(lt[f"{p}_y"], errors="coerce").values
                lcol = f"{p}_LLH"
                if lcol in ht.columns:
                    l = pd.to_numeric(ht[lcol], errors="coerce").values.astype(float)
                    llh[: min(N, len(l)), j] = np.nan_to_num(l[:N], nan=0.0)
            vis = (llh > 0.5) & np.all(np.isfinite(raw), axis=-1) & np.all(raw > 1, axis=-1)
            raw[~vis] = np.nan
            ud = np.full_like(raw, np.nan)
            ud[vis] = cams[k].undistort(raw[vis])
            raws.append(raw)
            uds.append(ud)
            sigs.append(np.where(llh < 0.85, 3.0, 1.0) * sig0[None, :])
        return raws, uds, sigs

    def set_offsets(self, offsets):
        self.offsets = np.asarray(offsets, int)
        self._build_arrays()

    def pidx(self, name):
        return self.points.index(name) if name in self.points else None

    def limit_mask(self):
        if self.limits is None:
            return np.ones(len(self.frames), bool)
        return (self.frames >= self.limits[0]) & (self.frames <= self.limits[1])

    def summary(self):
        lines = [f"{self.name}: cams {self.cam_ids}, offsets {self.offsets.tolist()}, limits {self.limits}, "
                 f"frames {self.frames[0]}..{self.frames[-1]} ({len(self.frames)}), res {self.W}x{self.H}",
                 f"  points present: {self.points}; missing: {self.missing_points}"]
        for k, c in enumerate(self.cam_ids):
            cnt = {p: int(self.vis[k, :, j].sum()) for j, p in enumerate(self.points)}
            lines.append(f"  cam{c}: " + " ".join(f"{p}={v}" for p, v in cnt.items()))
        return "\n".join(lines)
