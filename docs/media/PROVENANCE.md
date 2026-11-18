# Reconstruction preview

The README shows one excerpt from the camera-2 reconstruction of
`zorroBottomFlick2_20170903`. Video frames 33–62 (zero-based, end-exclusive 63)
correspond to absolute reconstruction frames 34–63. The MP4 retains all 30 frames
and plays at 9 frames/s, or 0.1× the 90-frames/s recording speed.

The existing CheetahHeadLab renderer projects the separate head and torso models
onto the camera image. The yellow arrow represents head-forward direction, not
measured eye gaze. The green box marks the annotated lure.

The source video includes the camera's fractional timing correction of
−0.2042194478673675 frames. The corrected overlay evaluates the model at the
absolute frame minus that correction. Its saved trajectory hash matches
`reference/zorroBottomFlick2_20170903/storedOptimRaw_trajectory.mat` in this release.

Source: `zorroBottomFlick2_20170903_recon_trajectory_pseudo_reproj_cam2.mp4`.
Its SHA-256 is
`3350a2d344310d2e5c5be2df689fbd3f5a7d07d4a788fed0a20137abd92ad90e`.
The reference trajectory SHA-256 is
`d663704eac0c3c173f44e1d20fe753881b40d55d38556958d00bfd084f725cc7`.

The excerpt was resized to 1280 × 720 and encoded as H.264 (CRF 19). The
720 × 405 GIF samples the same clip at 6 frames/s for a smaller looping preview.
The source framing, frame order, overlays and playback speed are preserved. No synthetic animal
imagery or interpolated video frames were added.

To regenerate the two files from the source video with FFmpeg:

```bash
ffmpeg -i zorroBottomFlick2_20170903_recon_trajectory_pseudo_reproj_cam2.mp4 \
  -vf 'trim=start_frame=33:end_frame=63,setpts=PTS-STARTPTS,scale=1280:720:flags=lanczos' \
  -an -c:v libx264 -preset slow -crf 19 -pix_fmt yuv420p \
  -movflags +faststart reconstruction.mp4
ffmpeg -i reconstruction.mp4 \
  -filter_complex '[0:v]fps=6,scale=720:405:flags=lanczos,split[a][b];[a]palettegen=max_colors=96:stats_mode=diff[p];[b][p]paletteuse=dither=bayer:bayer_scale=4' \
  -loop 0 reconstruction.gif
```

Original footage and calibration: [AcinoSet, African Robotics Unit](https://github.com/African-Robotics-Unit/AcinoSet).
Reconstruction overlays: CheetahHeadLab. Consult AcinoSet for the originating
recordings and their reuse terms.
