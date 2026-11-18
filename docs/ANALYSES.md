# Reproduce the manuscript analyses

Run these commands from the repository root. Reconstruction inputs are in
`data/inputs/`. Analysis outputs stay under the Git-ignored `outputs/` folder.
The original recordings are unnecessary for numerical analysis.

## Start from the included Monte Carlo results

```matlab
run_saved_analysis
```

This runs main and supplementary analyses from `data/analysis_cache/` without
reconstructing any poses or Monte Carlo draws. Then run the Python commands
in section 3 to render the panels and model comparison. Use
`run_saved_analysis(false)` to omit supplementary spatial and range-closure
analyses. The remaining sections describe the alternative workflow starting
from image annotations.

## 1. Reconstruct the pursuits

For a fresh all-pursuit reconstruction with the paper's uncertainty calculation:

```matlab
M = jsondecode(fileread('manifest.json'));
for k = 1:numel(M)
    run_reconstruction(M(k).dataset, ...
        'outputRoot','outputs/all18','nMC',500);
end
```

This is the expensive stage. For pitch statistics alone, nominal reconstruction
is sufficient: omit `'nMC',500` and skip datasets beginning with `zorro`.
If nominal outputs already exist, run `run_monte_carlo` on each output folder
instead of reconstructing it again. See [reproduction guide](REPRODUCING.md).

## 2. Calculate statistics

```matlab
run_pitch_analysis('outputs/all18')
run_yaw_analysis('outputs/all18')
```

With no argument, these commands read `outputs/<dataset>/`.
Both use `reference/manuscript_frames.csv`. Inputs are joined by absolute
frame ID, never by row number. Missing frames, failed uncertainty draws or
disagreement between nominal and Monte Carlo source reconstructions stop the
analysis. Prepared inputs are written to `outputs/analysis_inputs/`.

| Analysis | Cohort | Frames | Implementation |
|---|---|---:|---|
| Pitch ECDFs and centred distributions | 16 pursuits, excluding both close pass-by sequences | 1,590 | `PitchAnalysis/pitchTightnessAnalysis.m` |
| Pitch disturbance rejection and control gains | Same pitch cohort | 1,590 before finite-change windows | `PitchAnalysis/pitchControlModelAnalysis.m` |
| Yaw geometry and supported states | 11 frontal, 5 lateral, 2 close pass-by pursuits | 1,708 | `YawAnalysis/yawGeometryAndStates.m` |
| Absolute yaw-error ECDFs | 5 lateral and 2 close pass-by pursuits | 628 | `YawAnalysis/yawDistributionAnalysis.m` and ECDF renderers |
| Uncentred signed yaw distributions | 5 lateral pursuits | 510 | `YawAnalysis/plotYawSignedErrorDistributionLateralOnly.py` |
| Candidate aiming models | 5 lateral pursuits | 456 after the 100-ms velocity history | `YawAnalysis/compare_models.py` |

Pitch statistics use nominal orientations. Yaw ECDFs use the posterior-median
head and torso errors, wrapped to the shortest angular difference. Run-equal
summaries give each pursuit equal weight. Animal-equal ECDFs first average
runs within each cheetah, then give each cheetah equal weight.

For yaw tracking states, each joint draw uses the equal-range directional
projection `s = (h - rhat) dot u`. This distinguishes ahead from lagged head
orientation on either side of the target. A frame receives a supported state
when at least 95 percent of its 500 draws meet the joint angular and directional
criterion. Other frames remain unresolved. The main gate is 10 degrees.
The sensitivity figure uses 5 and 15 degrees.

## 3. Render panels and compare models

```bash
python -m pip install -r requirements.txt
python PitchAnalysis/plot_pitch.py
python YawAnalysis/plot_yaw.py
python YawAnalysis/compare_models.py
python tests/check_analysis.py
```

Pitch outputs include run- and animal-equal ECDFs, centred errors, continuous
disturbance rejection, control gains and a representative time series. Yaw
outputs include geometry/distance plots, context distributions, run- and
animal-equal ECDFs, signed errors and the supported-state/gate panels.
The export list is restricted to the [20 manuscript panels](FIGURE_COVERAGE.md).
Composite layouts are assembled separately. `python tests/check_figures.py`
checks the full export set after the supplementary analyses have run.
In MATLAB, `addpath('tests'); check_plot_scope` also verifies that direct calls
cannot enable the removed plots and that an empty spatial-example selection
still exports only the two examples shown in S15.

The model comparison evaluates current-target alignment, constant angular
lead, constant-velocity extrapolation and constant spatial lead. Each candidate
uses the same evaluation frames. Parameters are fitted on training pursuits
and evaluated on held-out pursuits, with an additional held-out-animal analysis.
The exported bar chart uses held-out-pursuit RMSE. No torso-referenced aiming
model is included.

`tests/check_analysis.py` checks cohort membership and compares numerical
outputs with the compact manuscript reference tables. Fresh reconstruction
may introduce small numerical differences. A mismatch is reported for inspection
rather than replaced by the stored result.

## Supplementary spatial lead and range closure

```matlab
run_yaw_analysis('outputs/all18',true)
```

This also runs the spatial-intersection and head/torso analyses for the five
lateral pursuits. Spatial lead is the forward head-ray/path intersection
measured along the reconstructed 3D target path. The endpoint is continued
along the terminal target direction, matching the updated supplementary
analysis. Invalid intersections remain missing. This spatial description is
distinct from the equal-range projection used for ahead/lagged classification.
Spatial calculations still use all five lateral pursuits. The plots comprise
the across-pursuit ratio panel and the two Monte Carlo examples used in S15.
The second example is displayed only before 0.64 s, matching the manuscript.

The two head/torso plots show range closure during yaw separation and a top-view
example (S16). Exploratory lagged-coupling results are retained as numerical tables
only. These diagnostics do not establish a
causal steering relationship. No turning pursuits enter the quantitative
manuscript cohort.

## Validation scope

The release analyses were run on the existing manuscript reconstructions and
500-draw ensembles. Pitch statistics, yaw frame errors, geometry, all three
state gates and the four candidate-model results matched the original tables
to numerical precision. Reconstruction itself was checked separately on the
examples described in [validation](VALIDATION.md). All 18 full Monte Carlo
ensembles were not regenerated during release preparation.
