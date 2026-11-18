# Yaw analysis

From the repository root, run `run_yaw_analysis` in MATLAB, then
`python YawAnalysis/plot_yaw.py`.

- `yawDistributionAnalysis.m`: paired posterior-median head/torso yaw errors and ECDFs.
- `yawGeometryAndStates.m` and `yawSampleGeometry.m`: target-travel geometry and Monte Carlo supported ahead/aligned/lagged states.
- `plot_yaw.py` and the accompanying renderers: geometry, distributions, ECDFs and alignment-gate sensitivity panels.
- `compare_models.py`: four candidate aiming rules with held-out-pursuit and held-out-animal validation.
- `yawSpatialLookAhead*`, `spatial*` and `plotSpatialLeadRatioAcrossPursuits.m`: forward head-ray/path intersections and reconstruction uncertainty.
- `headTorsoHandoff.m` and `plotHeadTorsoHandoff.m`: range closure, yaw separation and exploratory coupling diagnostics.

Geometry/state analyses use all 18 pursuits. The ECDF panel uses seven lateral
or close pass-by pursuits. Signed errors and candidate models use five lateral
pursuits. See the [analysis guide](../docs/ANALYSES.md) for frame counts and commands.
