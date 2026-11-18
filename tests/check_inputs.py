"""Check release checksums, camera indexing and frame coverage without MATLAB."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "CRBR/selfcalibration"))
from selfcalib.dataset import Dataset


def main():
    checksums = Path(__file__).with_name("SHA256SUMS")
    assert checksums.is_file(), f"Missing checksum inventory: {checksums}"
    for line in checksums.read_text().splitlines():
        digest, filename = line.split("  ", 1)
        actual = hashlib.sha256((ROOT / filename).read_bytes()).hexdigest()
        assert actual == digest, f"Checksum mismatch: {filename}"
    manifest = json.loads((ROOT / "manifest.json").read_text())
    assert len(manifest) == 18
    summary = []
    for entry in manifest:
        name = entry["dataset"]
        ds = Dataset(ROOT / "data/inputs" / name)
        assert ds.cam_ids == entry["camera_ids"], name
        assert ds.limits == (entry["first_frame"], entry["last_frame"]), name
        for labels, llh in zip(ds.label_tables, ds.llh_tables):
            # Some source visibility tables have unused trailing rows.
            assert len(llh) >= len(labels), f"Missing visibility rows: {name}"
        ref = loadmat(ROOT / "reference" / name / "storedOptimRaw_trajectory.mat")
        frames = ref["absFrames"].ravel().astype(int)
        assert len(frames) == entry["reference_frame_count"], name
        assert np.array_equal(frames, np.arange(ds.limits[0], ds.limits[1]+1)), name
        assert ref["xBarStore"].shape == (len(frames), 12), name
        assert np.isfinite(ref["xBarStore"]).all(), name
        for k, labels in enumerate(ds.label_tables):
            rows = frames + ds.offsets[k] - 1
            assert (rows >= 0).all() and (rows < len(labels)).all(), name
        if entry["original_input_included"]:
            raw = Dataset(ROOT / "data/original" / name)
            assert raw.cam_ids == ds.cam_ids, name
            assert not (ROOT / "data/original" / name / "selfcalib_report.json").exists()
        summary.append((name, ds.ncams, len(frames)))
    selected = ROOT / "reference/manuscript_frames.csv"
    if selected.exists():
        table = pd.read_csv(selected)
        assert len(table) == 1708 and table.pitch_included.sum() == 1590
        assert not table.duplicated(["dataset", "absolute_frame"]).any()
        for name, group in table.groupby("dataset"):
            ref = loadmat(ROOT / "reference" / name / "storedOptimRaw_trajectory.mat")
            assert np.isin(group.absolute_frame, ref["absFrames"].ravel()).all()
    for name, cameras, frames in summary:
        print(f"OK  {name}: {cameras} cameras, {frames} reconstruction frames")
    print(f"Verified inputs for {len(summary)} pursuits.")
    if selected.exists():
        print(f"Manuscript selection: {len(table):,} frames; {table.pitch_included.sum():,} for pitch.")
    print(f"Full reconstruction windows: {sum(x[2] for x in summary):,} frames before manuscript selection.")


if __name__ == "__main__":
    main()
