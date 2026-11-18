# Manuscript figure coverage

The standard workflow exports **20 numerical panels**, each used in the current
main manuscript or supplement. Run `run_saved_analysis` in MATLAB, then the
three Python rendering commands in the README. PDF and SVG versions are
alternative formats of the same panels. The scripts do not assemble or replace
the author's composite figures.

This inventory was checked against the current `paper/main.tex` and
`paper/supp.tex` in the development project (5 main figures and 17 supplementary
figures). Their hashes and the permitted export paths are recorded in
[`reference/manuscript_figures.json`](../reference/manuscript_figures.json).
Numbers refer to that manuscript version.

## Exported panels

Paths below are relative to `outputs/`, without the image extension.

| Manuscript panel | Content | Output |
|---|---|---|
| Figure 2c | Run-equal pitch ECDF | `pitch/tightness_run_equal/pitch_absolute_error_ecdf` |
| Figure 2d | Run-centred pitch errors | `pitch/tightness_run_equal/pitch_centered_error_distribution` |
| Figure 3b | Representative pitch coordination | `pitch/mechanism/pitch_representative_control_timeseries` |
| Figure 3c | Torso-to-head transmission gains | `pitch/mechanism/pitch_population_control_gains_run_level` |
| Figure 3d | Run-equal counter-rotation distribution | `pitch/mechanism/pitch_continuous_disturbance_rejection_run_equal` |
| Figure 4b | Seven-pursuit yaw ECDF | `yaw/yaw_absolute_error_ecdf_yaw_only` |
| Figure 4c | Five-lateral-pursuit signed yaw errors | `yaw/yaw_signed_error_distribution_lateral_only` |
| Figure 4d | Yaw lead versus pursuit geometry | `yaw/geometry/02_yaw_lead_vs_geometry` |
| Figure 4e | Yaw lead distributions by context | `yaw/geometry/01_yaw_lead_by_context` |
| Figure 4f | Monte Carlo-supported states, 10-degree gate | `yaw/yaw_mc_supported_modes_by_context` |
| Figure 5b | Candidate-model comparison | `yaw/models/candidate_model_comparison_lateral5` |
| Figure S10a | Animal-equal pitch ECDF | `pitch/tightness_animal_equal/pitch_absolute_error_ecdf` |
| Figure S10b | Animal-equal centred pitch errors | `pitch/tightness_animal_equal/pitch_centered_error_distribution` |
| Figure S13 | Animal-equal yaw ECDF | `yaw/yaw_absolute_error_ecdf_yaw_only_animal_equal` |
| Figure S14 | Supported states at 5- and 15-degree gates | `yaw/mode_threshold_sensitivity` |
| Figure S15a | Spatial lead divided by range | `yaw/spatial/supplementary_spatial_lead_ratio` |
| Figure S15b | Spatial lead and range with uncertainty | `yaw/spatial/MonteCarloResults/julesTopRun1_20171209_spatial_lookahead_monte_carlo_timeseries` |
| Figure S15c | Spatial lead and range, display before 0.64 s | `yaw/spatial/MonteCarloResults/julesTopRun11_20170829_spatial_lookahead_monte_carlo_timeseries` |
| Figure S16a | Closing speed versus head-torso separation | `yaw/head_torso/figure_6_closure_vs_decoupling` |
| Figure S16b | Representative top-view pursuit geometry | `yaw/head_torso/figure_7_top_view_decoupling_example` |

## Panels outside the current release plotting workflow

The release does **not yet regenerate every manuscript panel**. The remaining
items are listed here so that computed plots are not confused with hand-built
illustrations or missing exporters.

| Figure | Remaining material |
|---|---|
| 1 | Field photographs, model drawings, the CRBR flow diagram and pursuit-geometry illustrations. |
| 2a,b | Coordinate schematic and the representative image/3D/Monte Carlo pitch-error composite. The example's poses and draws are included, but its panel exporter and source footage are not. |
| 3a | Conceptual control diagram, including the proposed vestibular feedback loop. |
| 4a | Yaw-coordinate schematic. |
| 5a | Candidate aiming-rule illustration. |
| S1-S5 | Camera-consensus and optimisation illustrations, natural-pursuit photographs, synchronisation-viewer screenshot and epipolar-geometry illustration. |
| S6-S8 | Historical calibration and reconstruction diagnostics. Their diagnostic exporters and the required historical input snapshots are not bundled. |
| S9 | Current reconstruction-validation composite. Corrected inputs are included, but the composite exporter and original video frames are not. |
| S11-S12 | Lateral and close pass-by example composites. Saved poses and joint draws are included, but image overlays and example yaw-error panel exporters are not. |
| S17 | Qualitative turning-pursuit example. This separate sequence and its exporter are outside the release's 18-pursuit data bundle. |

S10 and S13-S16 have all their numerical panels covered by the export inventory.
S15c uses the manuscript's display limit before 0.64 s and a 0-40 m left axis.
This crop does not remove observations from the spatial statistics.

## Export checks

```bash
python tests/check_figures.py
python tests/check_analysis.py --include-spatial
```

The first check requires all 20 PNG panels and rejects unlisted PNG, PDF or SVG
exports under `outputs/pitch/` and `outputs/yaw/`. The second checks the numerical
results against the stored manuscript tables. Run the full supplementary
workflow before these checks. Files left by an earlier run are not deleted
automatically; use a fresh output folder or move those old plots elsewhere.

The restricted workflow omits the unused animal-equal disturbance plot, the
reference-adjusted diagnostic plot, the six nominal spatial time-series
exports, the three unused Monte Carlo spatial examples and their combined
plot, and the two duplicate yaw ECDF exports. Numerical tables, all five lateral
pursuits in the spatial calculations, and their uncertainty draws are retained.

## Direct plotting calls

The nominal spatial calculator exports tables only. Its old per-pursuit and
combined time-series renderers have been removed. The Monte Carlo spatial
renderer accepts only the two S15 example pursuits, including when called
directly or with `plotOnly=true`. Requesting the former combined plot or another
example stops with an explanatory error. An empty example list selects the two
manuscript examples, rather than all five pursuits.

The unused all-cohort yaw ECDF renderer and legacy pitch counter-rotation plot
have also been removed. Numerical tables remain available. Reconstruction
modules retain internal diagnostic visualisation helpers; these are disabled
by the documented release entry points and are not counted among the 20
manuscript analysis panels.
