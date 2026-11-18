# Optional reproducibility checks

Run these from the repository root after installing `requirements.txt`.
The analyses themselves do not depend on running these checks.

## Check the supplied files

```bash
python tests/check_inputs.py
```

This checks file fingerprints, camera numbering, annotation coverage, reference
trajectories and the exact manuscript frame selection. It does not reconstruct
poses or change data. `tests/SHA256SUMS` stores a SHA-256 fingerprint for each
supplied file, with paths relative to the repository root. A missing or changed
file fails the check. Deliberate edits also change a file's fingerprint, so this
check is intended for an unmodified download or clone.

The full reconstruction windows contain 1,900 frames. The analysis selects
1,708 of these, including 1,590 for pitch; the checker reports both counts.

## Check reproduced analyses

After `run_saved_analysis` and all three Python rendering commands:

```bash
python tests/check_analysis.py --include-spatial
python tests/check_figures.py
```

The first checks cohorts, numerical reference tables and Monte Carlo tracking
states. The second requires the 20 listed manuscript panels and rejects extra
analysis figure exports. Earlier output files are not deleted automatically.

To check the restrictions on individual plotting calls, run in MATLAB:

```matlab
addpath('tests')
check_plot_scope
```

## Check a fresh reconstruction

After reconstructing the default pursuit:

```matlab
addpath('tests')
check_reconstruction('outputs/julesBottomFlick2_20171209')
```

This compares the reconstructed trajectory against its saved numerical reference.
The tested examples and tolerances are described in
[validation](../docs/VALIDATION.md).
