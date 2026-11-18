# Reproduction guide

Run the commands below from the repository root. For a short introduction and
a reconstruction preview, see the [README](../README.md).

This package contains CRBR reconstruction, pitch and yaw analyses, and compact inputs for
*Cheetahs Combine Head Pitch Stabilisation with Prospective Yaw Orientation
during High-speed Pursuit*. It reconstructs separate head and torso models
and the lure from annotated multiview images, with optional joint Monte Carlo
reconstruction. Original videos are not needed for numerical reconstruction.

## Reproduce the manuscript panels

The shortest route uses the included nominal reconstructions and 500 joint
Monte Carlo draws. Follow the [README quick start](../README.md#reproduce-the-manuscript-panels)
to run `run_saved_analysis` and the three Python rendering commands. This exports
20 manuscript panels without repeating reconstruction. The sections below cover
optional reconstruction from image annotations.

## Reconstruct one pursuit

Open MATLAB in this folder and run:

```matlab
run_reconstruction
```

This reconstructs all 69 frames of `julesBottomFlick2_20171209` from the
released, self-calibrated annotations. It estimates scales, fits each frame,
reconstructs the lure, and applies the bidirectional trajectory estimator.
The saved production trajectory is used **only for comparison after fitting**.

Results are written to `outputs/julesBottomFlick2_20171209/`:

- `poses.csv`: absolute frame, head and torso rotations and origins.
- `angles_and_lure.csv`: lure positions and signed head/torso target errors.
- `storedOptimRaw_trajectory.mat`: the final MATLAB trajectory.
- `reference_comparison.csv`: framewise differences from the production pose.
- `reproduction_summary.json` and `reconstruction.log`: settings and checks.

The default example took about 3.3 minutes in the release validation. Runtime
depends on the machine and pursuit. To check agreement with the reference:

```matlab
addpath('tests')
check_reconstruction('outputs/julesBottomFlick2_20171209')
```

This example reproduced the reference exactly. A second, close-pursuit test
showed small pose differences and did not pass the same strict regression
tolerances. The measurements and interpretation are recorded in
`docs/VALIDATION.md`; exact agreement has not been established for all 18 runs.

Choose any of the 18 pursuits in `manifest.json`:

```matlab
run_reconstruction('julesTopRun1_20171209')
```

The runner refuses to overwrite an existing output folder. For another run
of the same example, choose a new root:

```matlab
run_reconstruction('julesBottomFlick2_20171209','outputRoot','outputs/repeat')
```

For every pursuit:

```matlab
M = jsondecode(fileread('manifest.json'));
for k = 1:numel(M)
    run_reconstruction(M(k).dataset,'outputRoot','outputs/all18');
end
```

## Requirements

Nominal reconstruction was tested with MATLAB R2026a on Linux. The retained
production code uses Optimization, Computer Vision, Symbolic Math, Signal
Processing and Image Processing Toolboxes. Behavioural analyses also use
Statistics and Machine Learning Toolbox. Parallel Computing Toolbox is
optional for parallel Monte Carlo execution. MATLAB's dependency report is
recorded in `docs/VALIDATION.md`. Earlier MATLAB releases have not been tested.
The supplementary spatial time-series exporter also uses Fontconfig (`fc-match`)
and MATLAB's Java runtime to position labels. Install Fontconfig for the tested
Linux workflow. Poppler's `pdftoppm` is recommended for PNG export; the renderer
falls back to MATLAB printing when it is unavailable.

Python runs self-calibration, input checks, figure rendering and yaw-model comparisons:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python tests/check_inputs.py
```

Use Python 3.10 or later. Tested dependency versions are recorded in
`docs/python_environment.txt`. `tests/check_inputs.py` verifies checksums, camera
ordering, annotation coverage, reference trajectories and manuscript frame IDs.

## Rerun spatiotemporal self-calibration

Corrected camera geometry, integer frame offsets and fractionally time-aligned
annotations are supplied for **all 18 pursuits**. These are the appropriate
inputs for the MATLAB quick start. Do not apply self-calibration a second time
to those corrected annotations.

To demonstrate the complete calibration stage, the package also includes the
original annotations, camera calibration and offsets for
`julesBottomFlick2_20171209`:

```bash
python run_selfcalibration.py
```

This runs the production temporal-block cross-validation and refinement of
camera timing and extrinsics. It writes a new corrected copy under
`outputs/selfcalibration/julesBottomFlick2_20171209/`, leaving the inputs intact.
Then reconstruct that newly corrected copy in MATLAB:

```matlab
run_reconstruction('julesBottomFlick2_20171209', ...
    'inputFolder','outputs/selfcalibration/julesBottomFlick2_20171209', ...
    'outputRoot','outputs/from_selfcalibration')
```

The other 17 pursuits are reproducible **from their supplied self-calibrated
inputs**. This minimal package does not include their pre-calibration inputs.
Each corrected dataset retains its production `selfcalib_report.json`, including
the estimated timing and camera corrections. Original source paths in these
reports are provenance only and are not used at runtime.

## Monte Carlo reconstruction

After nominal reconstruction:

```matlab
run_monte_carlo('outputs/julesBottomFlick2_20171209',500)
```

Or run both stages together in a new output folder:

```matlab
run_reconstruction('julesBottomFlick2_20171209', ...
    'outputRoot','outputs/with_uncertainty','nMC',500)
```

Each draw perturbs image annotations, model dimensions and lure annotations,
then refits the per-frame poses and the bidirectional trajectory. Camera
self-calibration remains fixed across draws, matching the production analysis.
The random seed is 1. All noise settings are explicit in
`CRBR/uncertainty/mcMainAnalysis.m`. Serial execution is the default.

For a quick execution check use `2` draws. Two draws do **not** provide usable
uncertainty intervals. The paper uses 500. Full Monte Carlo reconstruction is
much slower than the nominal example. Compact saved joint draws for all 18
pursuits are bundled in `data/analysis_cache/`; `run_saved_analysis` uses them
directly. Full optimizer diagnostics and intermediate fits are not bundled. To use MATLAB's parallel toolbox, pass `true` as the third argument
to `run_monte_carlo`, or `'useParallel',true` to `run_reconstruction`.

## What is included

| Folder/file | Purpose |
|---|---|
| `data/inputs/` | Corrected annotations, visibility scores, calibration, offsets and limits for 18 pursuits |
| `data/original/` | One original-input example for rerunning self-calibration |
| `reference/` | Small production pose/scale files for checking reconstructed results |
| `reference/manuscript_frames.csv` | Exact 1,708 analysed frame IDs, with the 1,590-frame pitch subset marked |
| `CRBR/` | Complete reconstruction pipeline described in the paper |
| `CRBR/selfcalibration/` | Camera timing and extrinsic refinement |
| `CRBR/reconstruction/` | Triangulation, camera selection, model scaling, pose fitting and bidirectional smoothing |
| `CRBR/lure/` | Lure triangulation only |
| `CRBR/uncertainty/` | Joint uncertainty reconstruction and angle calculations |
| `PitchAnalysis/` | Target-relative pitch precision and counter-rotation/control gains |
| `YawAnalysis/` | Yaw distributions, pursuit geometry, tracking states, candidate models and spatial lead |
| `docs/` | Input conventions, provenance and validation results |
| `tests/` | Optional input, output and reconstruction checks, including file checksums |

The full reconstruction windows contain 1,900 frames. They include boundary
and other frames outside the manuscript analysis cohort, since fitting and
smoothing precede behavioural frame selection. Use `manuscript_frames.csv`
when comparing outputs on the manuscript cohort. A raw reconstruction count
should not be reported as the number of analysed frames.

The [analysis guide](ANALYSES.md) runs the manuscript's pitch and yaw analyses
from reconstruction outputs. The scripts generate individual figure panels.
Manually assembled composite layouts, full source recordings, Blender projects
and intermediate reconstruction diagnostics are omitted. Compact Monte Carlo
analysis inputs are included. The README includes one short
reconstruction preview. Use a Git clone or GitHub source download to obtain the repository.
Generated outputs, virtual environments and caches are ignored by Git.
`tests/SHA256SUMS` verifies the supplied source files, inputs and reference data.

## Data provenance

Source recordings and camera calibration originate from
[AcinoSet](https://github.com/African-Robotics-Unit/AcinoSet). This package uses
the study's landmark annotations and refined camera geometry. See
`docs/DATA_AND_METHOD.md` and `docs/source_provenance.json` for file conventions
and source checksums. The README preview uses a short AcinoSet excerpt with reconstructed models
overlaid. Its source, timing and encoding are documented in
[clip provenance](media/PROVENANCE.md).
