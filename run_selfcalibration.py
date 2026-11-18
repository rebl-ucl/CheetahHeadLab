"""Run self-calibration on the released original-input example."""
from pathlib import Path
import argparse
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--out", type=Path,
                default=ROOT / "outputs/selfcalibration/julesBottomFlick2_20171209")
a = ap.parse_args()
subprocess.run([
    sys.executable, str(ROOT / "CRBR/selfcalibration/scripts/run_selfcalib.py"),
    str(ROOT / "data/original/julesBottomFlick2_20171209"), str(a.out),
    "--cv", "--fps", "90",
], check=True)
