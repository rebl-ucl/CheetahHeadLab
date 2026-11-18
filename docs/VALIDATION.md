# Release validation

Prepared on 17 November 2026. Tests use the staged standalone package and a
clean MATLAB search path (`restoredefaultpath`). Core project folders are not
on the test path. Inputs are copied to a fresh output directory before each
reconstruction, so no saved pose or triangulation is used as an input.

## Input inventory

`python tests/check_inputs.py` passed for all 18 pursuits, including camera order,
annotation/visibility coverage, reconstruction limits and reference frame
coverage. The reconstruction windows contain 1,900 frames. The independently
listed manuscript cohort contains 1,708 unique dataset/frame pairs, with
1,590 marked for pitch analysis.

## Nominal reconstruction

The five-camera `julesBottomFlick2_20171209` example reconstructed all 69
frames from corrected image annotations in 198.49 seconds. Frame coverage
matched the production reference. Maximum differences in head/torso forward
direction and origin position were all zero on this machine.

The six-camera close pursuit `zorroBottomFlick2_20170829` reconstructed all
35 frames in 112.32 seconds. Triangulated landmarks, camera masks and model
scales matched the saved production intermediates exactly. The subsequent
pose fit was not bit-identical: maximum head-axis difference was 0.00089
degrees and maximum torso-axis difference was 0.1232 degrees. Head-origin
difference was at most 0.0016 mm and torso-origin difference was 15.47 mm.
The strict regression test therefore **did not pass for this second run**.

Comparing the fitted anatomical torso landmarks gives much smaller maximum
differences: 0.553 mm at the neck and 0.186 mm at the spine. The torso model
has two landmarks along one axis, so rotation about that axis is weakly
constrained while its off-axis origin can move. This is consistent with the
observed parameter differences, but the exact numerical cause has not been
established. The discrepancy is recorded rather than concealed by relaxing
the regression thresholds. See `zorro_landmark_comparison.json` and the
nominal comparison summaries supplied here.

## Monte Carlo execution check

Two independent joint draws were run across all 69 frames of the default
example with seed 1 and serial execution. Both succeeded and all pose
parameters were finite. The nominal reconstruction rerun through the Monte
Carlo code path agreed with the saved nominal trajectory to less than
1e-6 degrees. This is an execution check, not a substitute for the 500-draw
uncertainty calculation used by the manuscript.

## Self-calibration

The original-input example was rerun with four-fold temporal-block
cross-validation and the final full-data fit. It completed in 465 seconds.
As in the production report, it selected joint timing and extrinsic refinement.
Integer camera offsets matched exactly. The largest timing-correction
difference was 0.00384 frames, and the largest individual corrected image
coordinate difference was 0.440 pixels. Camera intrinsic parameters and
distortion coefficients were unchanged. Thus this rerun reproduced the
selected method and a close calibration solution, not bit-identical
corrected inputs. The numerical comparison is in
`selfcalibration_validation.json`.

## Environment and dependencies

MATLAB: R2026a, version 26.1.0.3203278, Linux. Python package versions are in
`python_environment.txt`.

MATLAB's static dependency analysis identified MATLAB, Optimization, Signal
Processing, Symbolic Math, Image Processing, Parallel Computing and Computer
Vision Toolboxes. Parallel Computing Toolbox is only needed when parallel
Monte Carlo execution is requested. The wrapper uses serial execution by
default. The static code dependency list was used to select the MATLAB
source files retained from the production pipeline. The behavioural analyses
also require Statistics and Machine Learning Toolbox.

`tests/check_reconstruction.m` compares every reconstructed frame with its
reference and checks finite poses. Its tolerances (0.01 degrees for forward
axes and 0.1 mm for origins) are numerical regression criteria, not claims
about biological accuracy or measurement uncertainty.

The full 500-draw ensembles for all 18 pursuits were not rerun as part of
release validation. They can be generated using the supplied code and inputs.

## Pitch and yaw analyses

The reorganised release was tested with a clean MATLAB search path against
the existing manuscript reconstructions and 500-draw ensembles. No core code
folders were added to the path. The released input adapter reproduced the
stored kinematic quantities with a maximum absolute difference of 8.89e-16.

The following outputs agreed with the manuscript tables within 1e-8 numerical
tolerance: pitch ECDFs and centred distributions, 240 pitch run/delay control
fits, all 1,708 paired yaw errors, the 18 pursuit geometry summaries, and all
5,124 frame/gate probability rows at 5, 10 and 15 degrees. Four candidate yaw
models reproduced the held-out-pursuit and held-out-animal results across
456 evaluation frames from the five lateral pursuits.

The spatial-lead analysis and its 500-draw uncertainty propagation completed
for all five lateral pursuits. Their summary tables matched the updated
supplementary results, including nominal median d/R values of 0.24 to 0.56.
The range-closure and exploratory head/torso analyses also completed.

The Python renderers were executed and representative output panels were
visually inspected. `tests/check_analysis.py` provides regression checks
against compact reference tables. These checks validate the extracted analysis
code on the existing reconstructions, independently of the nominal and
two-draw reconstruction checks above.

## Manuscript-only figure exports

The full saved-data workflow was rerun in an isolated copy after restricting
plot exports to the current manuscript inventory. It produced 20 PNG panels,
with PDF/SVG alternatives and no unlisted figure exports. All existing analysis
checks passed, including spatial lead. Numeric values in all 62 CSV tables shared
with the previous complete run agreed within 1e-9 absolute and relative tolerance.
The export checker also rejected intentionally added extra plots and missing panels.

This validates the export selection and unchanged numerical results. It does
not fill the missing example and validation exporters listed in
[figure coverage](FIGURE_COVERAGE.md).

Direct-call checks now also cover the removed export paths. Requests for
nominal-only spatial plots, a combined Monte Carlo spatial plot, a spatial
example absent from S15, or the old all-cohort yaw ECDF plot were all rejected
before loading data. A plot-only call with an empty example list exported
exactly the two S15 examples. The full workflow was rerun after these changes:
20 panels, no extra exports, and all 62 shared numerical CSV tables unchanged
within the tolerance above.
