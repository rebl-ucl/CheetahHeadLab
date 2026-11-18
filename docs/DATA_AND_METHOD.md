# Data and method conventions

## Inputs

- `camN_labels_wide.csv` contains image coordinates in pixels. Camera N is
  matched to the corresponding position in the sorted camera-ID list.
- `camN_llh_wide.csv` contains visibility/confidence scores. Some source tables
  have extra trailing rows beyond the coordinate table; these rows are unused.
- `AcinoSet_data/<n>_cam_scene_sba.json` contains the camera intrinsics,
  distortion coefficients and extrinsics.
- `AcinoSet_data/offset.xlsx` stores integer camera offsets by camera-list
  position. For absolute frame f, the MATLAB row is f + offset. Python uses
  zero-based row f + offset - 1. Preserve the full annotation tables and their
  row ordering even when only part of the recording is reconstructed.
- `limits.xlsx` gives the inclusive reconstruction interval. The frame rate is
  taken explicitly from `manifest.json`; no video-based inference is needed.
- `selfcalib_report.json` records the production calibration corrections.
  Fractional corrections are already applied to the supplied image labels.

Only the seven model landmarks are used for head/torso fitting: left/right
eyes, nose, left/right ears, spine and neck base. Other annotated landmarks
remain in the source tables because self-calibration can use them. Lure
annotations feed a separate triangulation step.

## CRBR reconstruction stages

CRBR names the complete pipeline, including spatiotemporal self-calibration
and multiview triangulation. Its implementation is grouped under `CRBR/`.

1. **Spatiotemporal self-calibration.** The Python stage fits timing and camera
   extrinsics and exports camera geometry, integer offsets and fractionally
   time-aligned annotations. The original-input example includes temporal-block
   cross-validation using the same production script.
2. **Multiview triangulation.** `multiviewRecon` undistorts image coordinates
   and applies the leave-one-camera-out camera-trust procedure.
3. **Head and torso model scaling.** `estimateScales` obtains the two model
   scales from accepted multiview landmarks.
4. **Nonlinear least-squares optimisation.** `reconPerFrame` fits the head and
   torso model poses using robust, consensus-weighted image residuals.
5. **Bidirectional constant-velocity smoothing.** `reconTrajectoryBatch`
   repeatedly refits interior frames using the observations and robust
   neighbouring-frame rotational and translational priors. Rotations are
   regularised on SO(3). `lureReconTriangulation` supplies the target positions.
6. **Joint Monte Carlo reconstruction.** Optional draws repeat pose and lure
   reconstruction with perturbed labels and model dimensions.

The wrapper retains production numerical defaults: Huber threshold 3,
2-degree rotational midpoint sigma, 0.015-m translational midpoint sigma,
temporal Huber threshold 3, and up to six trajectory passes. The full source
contains the remaining tolerances and stopping rules.

## Output conventions

`poses.csv` contains 12 pose parameters. Head/torso origins are in world metres.
Angles are radians with the production rotation convention
`R = Rz(rz) * Rx(rx) * Ry(ry)`, using the rotation matrices in
`CRBR/uncertainty/mcMainAngles.m`. Do not substitute a library's default
Euler convention. The forward vector is `R * [1; 0; 0]`.

`angles_and_lure.csv` reports degrees. Alpha/beta are head/torso pitch errors
relative to the target. Gamma/delta are wrapped azimuth differences relative
to the current target bearing. Their signs are coordinate signs, **not** the
travel-directed ahead/lag classification in the manuscript. Missing target
positions produce missing target-error angles while retaining the pose.

Saved production poses are kept outside the input folders. The wrapper never
uses them to initialise or fit a reconstruction. It reports differences in
forward-axis direction and origin position, because different Euler triples
can describe equivalent orientations.

`reference/manuscript_frames.csv` records the final frame selection from the
existing analysis outputs and quality windows. Both close pass-by pursuits
are excluded from its pitch subset. `prepare_analysis_inputs` applies these
exact IDs after reconstruction, verifies pose/lure/Monte Carlo agreement and
calculates world-frame and target-relative angles. The fixed frame selection
is specific to the supplied pursuits. New recordings require their own quality
assessment. Yaw-state classification is implemented in
`YawAnalysis/yawGeometryAndStates.m`.

## Source snapshot

The numerical MATLAB functions and Python self-calibration modules originate
from the working pipeline. In this release, they are organised under `CRBR/`
with updated relative paths and diagnostic labels. The
self-calibration CLI adds an explicit frame-rate option, rejects existing
output folders, and disables video copying. These changes do not alter the
fitting method. `docs/source_provenance.json` records hashes of both versions
where a release adaptation exists. Its commit ID is contextual because the
source working tree also contains uncommitted changes.

This package does not assign a new licence to the originating AcinoSet data.
Users should consult the originating dataset's terms when reusing those data.
