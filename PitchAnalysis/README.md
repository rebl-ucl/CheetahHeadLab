# Pitch analysis

From the repository root, run `run_pitch_analysis` in MATLAB, then
`python PitchAnalysis/plot_pitch.py`.

- `pitchTightnessAnalysis.m`: absolute-error ECDFs and run-centred error distributions, with run- and animal-equal summaries.
- `pitchControlModelAnalysis.m`: finite-change disturbance rejection, counter-rotation and control gains.
- `plot_pitch.py`: manuscript panel styling from the exported numerical tables.

The cohort contains 16 pursuits and 1,590 frames. Both close pass-by pursuits
are excluded. Inputs come from CRBR outputs, selected using the fixed manuscript
frame list. See the [analysis guide](../docs/ANALYSES.md) for the full workflow.
