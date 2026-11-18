# Saved manuscript reconstruction ensembles

These compact inputs contain the nominal head/torso trajectories, lure paths,
and 500 successful joint Monte Carlo pose/lure draws for each of 18 pursuits.
Monte Carlo arrays retain the 1,708 manuscript frames, with absolute frame IDs,
posterior-median angles and source reconstruction metadata. The nominal files
retain their original reconstruction coverage. Both close pass-by pursuits
are excluded automatically from pitch analysis.

Run `run_saved_analysis` from the repository root in MATLAB. This recomputes
statistics from these inputs without refitting the Monte Carlo reconstructions.
Outputs go to `outputs/`; these saved inputs are read only.

The ensembles originate from the manuscript reconstruction outputs in
CheetahHeadLab_core, extracted by prepare_analysis_inputs and validated against
manuscript tables. Original paths inside MATLAB metadata record provenance
only and are not required at runtime. Per-file hashes are in manifest.json
and the repository SHA256SUMS. Optimizer diagnostics and intermediate fits
are omitted; all joint poses and lure positions needed by the analyses remain.
