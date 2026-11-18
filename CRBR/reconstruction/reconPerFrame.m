% reconPerFrame.m
%
% CRBR per-frame pose reconstruction. The head model uses left/right eyes,
% nose and left/right ears. The torso uses spine and neck-base landmarks.
% The kinematic model (cheetah3D) and reprojection error (cheetahError)
% fit these landmarks jointly.
%
% This stage:
%   1) Uses the multiview triangulation results from multiviewRecon.m
%      (triangulatedPointsRaw_multiview.mat -> TRI.pointsRaw, which now also
%      contains l_ear/r_ear) as the SEED for EVERY reconstruction frame.
%   2) Disables the temporal smoothness term in cheetahError for all
%      frames (useTemporal = false). Each frame is fit purely to its own
%      multi-camera reprojection evidence. A robust loss is available as an
%      option: robustOpt is passed to the lsqnonlin solves below and selects
%      the per-observation outlier handling inside cheetahError. It
%      defaults to 'none' (plain least squares, unweighted default) and
%      can be switched to 'huber'/'tukey'/'hard'.
%
% Missing triangulated landmarks at a given frame (NaN entries in
% pointsRaw) are filled with the previous frame's kinematic-model points
% before the seed-refinement step. The main optimiser sees only
% reprojection residual, so the fallback acts only as an initial guess.
%
% Heading "addition" (0 or 180 deg) is derived from the multiview body
% axis (spine -> neck_base) per frame, not from AS.x_0. AS is loaded
% only as a frame-1 fallback if TRI is missing.
%
% Outputs are written to:
%   <parent_folder>/reconResults_CRBR_perFrame/
%   <parent_folder>/storedOptimRaw.mat
%   <parent_folder>/storedOptimMed.mat
%   <parent_folder>/storedOptimLPF.mat
%   <parent_folder>/cheetahReprojErr.mat

function reconPerFrame(parent_folder, camIDs, opts)
% reconPerFrame  CRBR per-frame head and torso pose estimation.
%
%   reconPerFrame(parent_folder, camIDs, opts)
%
% CRBR per-frame reconstruction: seeds every frame from the multiview
% triangulation (triangulatedPointsRaw_multiview.mat, including ears) and fits
% the 12-DOF head+torso pose to each frame's reprojection evidence with no
% temporal smoothing. Writes storedOptimRaw.mat / storedOptimMed.mat /
% storedOptimLPF.mat (+ diagnostics) to <parent_folder>.
%
% opts fields (all optional; defaults are headless + fast for batch):
%   showFigures(false)  saveFigures(true)  createVideo(false)  createOrthoVideo(false)
%   videoSource("raw")  diagnose_frame_1(false)  show_frame1_init_vs_sol(false)
%   optimSuffix("")  tag appended to the saved storedOptim*/cheetahReprojErr
%                    filenames (keep "" for the default reconstruction)
%
% Called with no arguments it runs a default dataset (standalone debugging).

if nargin < 1 || isempty(parent_folder)
    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    parent_folder = fullfile(repoRoot, 'outputs', 'phantomTopRun1_20170903');
end
if nargin < 2 || isempty(camIDs)
    camIDs = [1 2 3 4 5 6];
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

ncams  = numel(camIDs);

% -------------------- Optimisation / Plot / Filter options --------------------
% Figure / video show/save convention (shared across the pipeline):
%   showFigures - pop up diagnostic figure windows
%   saveFigures - write the diagnostic report PNGs to reconResults_CRBR_perFrame/
%   showVideos  - render the reconstruction video live while writing
%   saveVideos  - write the reconstruction video (createVideo) to .avi
showFigures             = getOpt(opts, 'showFigures',             false);
saveFigures             = getOpt(opts, 'saveFigures',             true);
showVideos              = getOpt(opts, 'showVideos',              false);
saveVideos              = getOpt(opts, 'saveVideos',              false);
createVideo             = getOpt(opts, 'createVideo',             saveVideos);   % primary 3D recon video
createOrthoVideo        = getOpt(opts, 'createOrthoVideo',        false);        % secondary ortho video (opt-in)
videoSource             = getOpt(opts, 'videoSource',             "raw");   % "raw", "med", or "lpf"
diagnose_frame_1        = getOpt(opts, 'diagnose_frame_1',        false);   % frame-1 diagnostic figures
show_frame1_init_vs_sol = getOpt(opts, 'show_frame1_init_vs_sol', false);   % frame-1 init-vs-solution overlay
% Build the diagnostic figures whenever they will be shown or saved; each
% figure's on-screen visibility tracks showFigures, the PNG write tracks
% saveFigures (so unattended batch runs save without popping windows).
doFigs                  = showFigures || saveFigures;
% Optional tag appended to the saved storedOptim*/cheetahReprojErr filenames.
% Default "" keeps the canonical storedOptim{Raw,Med,LPF}.mat names.
optimSuffix             = string(getOpt(opts, 'optimSuffix', ""));

resultsSubdir = "reconResults_CRBR_perFrame";

if ~exist(fullfile(parent_folder, resultsSubdir), 'dir')
    mkdir(fullfile(parent_folder, resultsSubdir));
end

% -------------------- Load geometry scales --------------------
scaleFile = fullfile(parent_folder, "cheetah_scales.mat");

if exist(scaleFile, 'file') == 2
    Sscale = load(scaleFile, "scaleHead", "scaleTorsoLength");
    scaleHead        = Sscale.scaleHead;
    scaleTorsoLength = Sscale.scaleTorsoLength;
    fprintf("Loaded scales from cheetah_scales.mat\n");
else
    warning("Scale file not found. Using default scales = 1.");
    scaleHead        = 1;
    scaleTorsoLength = 1;
end

fprintf("scaleHead        = %.4f\n", scaleHead);
fprintf("scaleTorsoLength = %.4f\n", scaleTorsoLength);

% -------------------- Load manual labels + LLH CSVs --------------------
for ii = 1:ncams
    cam = camIDs(ii);
    eval(sprintf('cam%dlabels = readtable(fullfile(parent_folder, "cam%d_labels_wide.csv"));', cam, cam));
    eval(sprintf('cam%dLLH    = readtable(fullfile(parent_folder, "cam%d_llh_wide.csv"));',    cam, cam));
end

for ii = 1:ncams
    i = camIDs(ii);
    labels = eval(sprintf('cam%dlabels', i));
    llh    = eval(sprintf('cam%dLLH', i));
    fprintf("height(cam%dlabels)=%d, height(cam%dLLH)=%d\n", i, height(labels), i, height(llh));
end

% -------------------- Body-part names (ears included for CRBR) --------------------
indexes = struct();
indexes.l_eye     = [];
indexes.r_eye     = [];
indexes.nose      = [];
indexes.l_ear     = [];
indexes.r_ear     = [];
indexes.spine     = [];
indexes.neck_base = [];

idxs = struct2table(indexes);
varnames = idxs.Properties.VariableNames;

varsUsed = {'l_eye', 'r_eye', 'nose', 'l_ear', 'r_ear', 'spine', 'neck_base'};

% -------------------- Offsets --------------------
osfile = dir(fullfile(parent_folder, "AcinoSet_data", "*offset*"));
if ~isempty(osfile)
    offsets = readtable(fullfile(parent_folder, "AcinoSet_data", osfile.name), ...
        VariableNamingRule="preserve");
else
    offsets = table(camIDs(:), zeros(ncams,1), 'VariableNames', {'cam','offset'});
end

vals = offsets.offset;
for k = 1:numel(vals)
    fprintf('offset(%d) = %g\n', k, vals(k));
end

% -------------------- Load camera parameters --------------------
jsonCandidates = dir(fullfile(parent_folder, "AcinoSet_data", "*cam_scene_sba.json"));
assert(~isempty(jsonCandidates), "No *cam_scene_sba.json found.");

wantedName = sprintf("%d_cam_scene_sba.json", ncams);
match = find(strcmp({jsonCandidates.name}, wantedName), 1);
assert(~isempty(match), "Could not find %s in AcinoSet_data.", wantedName);

extrinsicsFile = fileread(fullfile(parent_folder, "AcinoSet_data", jsonCandidates(match).name));
extrinsics     = jsondecode(extrinsicsFile);

poses           = {};
intrinsicParams = {};

for ii = 1:ncams
    c     = extrinsics.cameras(ii);
    R_cw  = c.r;
    t_cw  = c.t';
    R_wc  = R_cw.';
    C_w   = -(R_cw.' * t_cw.').';
    T_wc  = rigidtform3d(R_wc, C_w);
    poses{ii} = T_wc;

    imageSize = [extrinsics.camera_resolution(2) extrinsics.camera_resolution(1)];
    camIntrin = cameraIntrinsicsFromOpenCV(c.k, c.d, imageSize);
    intrinsicParams{ii} = camIntrin;
end

if ncams == 6
    plotCameraConfig6(poses, 0.1);
else
    plotCameraConfigLT6(poses, 0.1, camIDs);
end
hold on;

% -------------------- Load AS (used only as fallback for frame 1 heading) --------------------
hasAS = false;
getAS = dir(fullfile(parent_folder, "AcinoSet_data", "*positions3DwithX.mat"));
if ~isempty(getAS)
    AS = load(fullfile(parent_folder, "AcinoSet_data", getAS(1).name));
    if isfield(AS, "x") && isfield(AS, "positions")
        states = AS.x;
        hasAS  = true;
        fprintf("AS data loaded (fallback only). height(AS.x)=%d, size(AS.positions)=%s\n", ...
            height(AS.x), mat2str(size(AS.positions)));
    end
end
if ~hasAS
    warning("AS data not available. Frame-1 heading will use multiview geometry only.");
    states = [];
end

% -------------------- Reconstruction frame range from Blender LABELS --------------------
refCam    = camIDs(1);
refLabels = eval(sprintf('cam%dlabels', refCam));
if ismember("frame", refLabels.Properties.VariableNames)
    refLabels = removevars(refLabels, "frame");
end

limitfile = dir(fullfile(parent_folder, "*limit*"));
if ~isempty(limitfile)
    lims = readtable(fullfile(parent_folder, limitfile(1).name), ...
        VariableNamingRule="preserve");
    startI = lims{1,1};
    endI   = lims{1,2};
    fprintf("Reconstruction range from limit file: abs frame %d -> %d\n", startI, endI);
else
    startI = 1;
    endI   = height(refLabels);
    fprintf("No limit file. Defaulting to all available frames: abs frame %d -> %d\n", startI, endI);
end

% Honour the declared acquisition range.  Older code forced startI >= 10,
% which silently shifted/truncated datasets whose limits begin at frames 1--9
% while triangulation, lure reconstruction, and downstream analyses used the
% actual limit-file start.
startI = max(1, startI);
endI   = min(height(refLabels), endI);
assert(startI <= endI, "Invalid reconstruction range.");

nframes = endI - startI + 1;
absFrames = (startI:endI).';
fprintf("Final range: %d -> %d (nframes=%d)\n", startI, endI, nframes);

% -------------------- Load + undistort Blender labels --------------------
undistortedLabels = struct();
likelihoods       = struct();

for ii = 1:ncams
    cam       = camIDs(ii);
    camLabels = eval(sprintf('cam%dlabels', cam));
    LLHtbl    = eval(sprintf('cam%dLLH', cam));

    if ismember("frame", camLabels.Properties.VariableNames)
        camLabels = removevars(camLabels, "frame");
    end
    if ismember("frame", LLHtbl.Properties.VariableNames)
        LLHtbl = removevars(LLHtbl, "frame");
    end

    assert(ismember("cam", offsets.Properties.VariableNames), ...
        "offsets table must contain a cam column.");
    row = find(offsets.cam == cam, 1);
    assert(~isempty(row), "Camera %d not found in offsets.", cam);
    offset = offsets.offset(row);

    intrin   = intrinsicParams{ii};
    rowStart = startI + offset;
    rowEnd   = endI   + offset;
    assert(rowStart >= 1, "cam%d: rowStart < 1", cam);
    assert(rowEnd <= height(camLabels), "cam%d: rowEnd exceeds height(camLabels)", cam);
    assert(rowEnd <= height(LLHtbl),    "cam%d: rowEnd exceeds height(LLHtbl)",    cam);

    camLabelsCrop = camLabels(rowStart:rowEnd, :);
    LLHcrop       = LLHtbl(rowStart:rowEnd, :);

    [camUDx, camUDy] = undistortBlenderLabels(camLabelsCrop, intrin, nframes, varnames);

    camUD = camUDx;
    for k = 1:width(camUDx)
        xVals = camUDx{:, k};
        yVals = camUDy{:, k};
        camUD.(camUDx.Properties.VariableNames{k}) = ...
            arrayfun(@(x, y) [x, y], xVals, yVals, 'UniformOutput', false);
    end

    camname = sprintf('cam%d', cam);
    undistortedLabels.(camname) = camUD;

    wantedLLH = {'l_eye_LLH','r_eye_LLH','nose_LLH','l_ear_LLH','r_ear_LLH','spine_LLH','neck_base_LLH'};
    LLHcrop   = LLHcrop(:, wantedLLH);
    LLHcrop.Properties.VariableNames = varsUsed;

    likelihoods.(camname) = LLHcrop;
end

% -------------------- Load multiview triangulation results (PRIMARY seed source) --------------------
triFile = fullfile(parent_folder, "triangulatedPointsRaw_multiview.mat");
assert(isfile(triFile), ...
    "Missing %s. Run multiviewRecon.m first.", triFile);

TRI = load(triFile);
assert(isfield(TRI, 'pointsRaw') && isfield(TRI, 'absFrames'), ...
    "%s must contain pointsRaw and absFrames.", triFile);

% Build a map from absolute frame number -> row in TRI.pointsRaw
% TRI.pointsRaw.(landmark) is [nTRI x 3]; TRI.absFrames is [nTRI x 1].
% Ears (l_ear, r_ear) must be present for the CRBR reconstruction.
landmarkNames = {'l_eye','r_eye','nose','l_ear','r_ear','spine','neck_base'};
for k = 1:numel(landmarkNames)
    assert(isfield(TRI.pointsRaw, landmarkNames{k}), ...
        "pointsRaw missing field '%s'. Re-run multiviewRecon.m so ears are triangulated.", landmarkNames{k});
end

triAbs   = TRI.absFrames(:);
triRowOf = containers.Map('KeyType','double','ValueType','double');
for r = 1:numel(triAbs)
    triRowOf(double(triAbs(r))) = r;
end

fprintf("Loaded multiview triangulation: %d frames (abs %d..%d)\n", ...
    numel(triAbs), min(triAbs), max(triAbs));

% -------------------- Per-frame camera-trust mask (from multiview LOO pass) --------------------
% multiviewRecon.m flagged cameras carrying a systematic bias per frame via
% leave-one-out reprojection (TRI.camKeptMask, rows aligned to TRI.absFrames). Map
% it onto this run's fit frames; a false entry HARD-EXCLUDES that camera from the
% frame's cheetahError fit (same as LLH=0). Default = keep all cameras, so a
% triangulation .mat without the field reproduces the previous behaviour.
camMaskStore = true(nframes, ncams);
camLOOStore  = nan(nframes, ncams);   % per-frame LOO residual carried over for the log
if isfield(TRI, 'camKeptMask') && isfield(TRI, 'camIDs')
    triCamIDs = double(TRI.camIDs(:)).';
    [~, colOfLocal] = ismember(camIDs, triCamIDs);   % this run's cams -> TRI columns
    hasLOO = isfield(TRI, 'camLOOResid');
    for i = 1:nframes
        a = startI + (i - 1);
        if isKey(triRowOf, double(a))
            row = triRowOf(double(a));
            for c = 1:ncams
                if colOfLocal(c) > 0
                    camMaskStore(i,c) = logical(TRI.camKeptMask(row, colOfLocal(c)));
                    if hasLOO
                        camLOOStore(i,c) = TRI.camLOOResid(row, colOfLocal(c));
                    end
                end
            end
        end
    end
    fprintf("Per-frame camera-trust mask loaded: %d camera-frames rejected (of %d).\n", ...
        sum(~camMaskStore(:)), numel(camMaskStore));
else
    fprintf("No camKeptMask in triangulation file; using all cameras in the fit.\n");
end

% (Local functions getMVseed, headingFromMVorPrev, buildParams0 defined at end of file.)

% -------------------- Kinematic chain model --------------------
syms theta_x theta_y theta_z theta_x_2 theta_y_2 theta_z_2
syms x y z x2 y2 z2

Rx = @(theta_x) [1 0 0; 0 cos(theta_x) sin(theta_x); 0 -sin(theta_x) cos(theta_x)];
Ry = @(theta_y) [cos(theta_y) 0 sin(theta_y); 0 1 0; -sin(theta_y) 0 cos(theta_y)];
Rz = @(theta_z) [cos(theta_z) sin(theta_z) 0; -sin(theta_z) cos(theta_z) 0; 0 0 1];

params = [theta_x theta_y theta_z x y z theta_x_2 theta_y_2 theta_z_2 x2 y2 z2];

% -------------------- Frame-1 seed from multiview triangulation --------------------
[seedMV1, validMask1] = getMVseed(startI, triRowOf, TRI, landmarkNames);
nValid1 = sum(validMask1);
fprintf("Frame 1 multiview seed: %d / %d landmarks valid.\n", nValid1, numel(landmarkNames));

assert(nValid1 >= 2, ...
    "Multiview triangulation has <2 valid landmarks at frame %d; cannot seed frame 1.", startI);

% Frame-1 prev_pts fallback: use AS landmarks at the corresponding AS
% index if AS is available (so a missing MV landmark at frame 1 falls
% back to AS rather than to [0 0 0]). Otherwise use zeros and rely on
% multiview coverage to be sufficient. AS has no ears, so ear fallbacks
% remain [0 0 0] and the ears must come from the multiview seed.
prev_pts_for_frame1 = struct( ...
    'l_eye',     [0 0 0], 'r_eye',     [0 0 0], 'nose',      [0 0 0], ...
    'l_ear',     [0 0 0], 'r_ear',     [0 0 0], ...
    'spine',     [0 0 0], 'neck_base', [0 0 0]);
% AS (AcinoSet) frame-1 seed kept for the init/solution/AS 3D overlay below:
%   asSeedIdx    - AS.positions slice index for this frame (NaN if no AS); the
%                  FULL AcinoSet landmark cloud is plotted straight from
%                  AS.positions(:,:,asSeedIdx), matching the "AS seed" figures
%                  in the CRBR pose model.
%   asPts_frame1 - the head/neck AS landmarks as a labelled struct (finite only),
%                  used only to draw an AS orientation skeleton.
asSeedIdx    = NaN;
asPts_frame1 = struct();
if hasAS
    AS_abs_start_path = fullfile(parent_folder, "AcinoSet_data/", "FTE_start_frame.xlsx");
    if isfile(AS_abs_start_path)
        T_AS         = readtable(AS_abs_start_path);
        AS_abs_start = T_AS{1,2};
        asIdx_f1     = max(1, min(size(AS.positions,3), startI - AS_abs_start + 1));
        asSeedIdx    = asIdx_f1;
        asRows       = [1 2 3 4 5];  % l_eye, r_eye, nose, neck_base, spine in AS
        asFieldOrder = {'l_eye','r_eye','nose','neck_base','spine'};
        for k = 1:numel(asFieldOrder)
            v = squeeze(AS.positions(asRows(k), :, asIdx_f1))';
            v = v(:).';
            if all(isfinite(v))
                prev_pts_for_frame1.(asFieldOrder{k}) = v;
                asPts_frame1.(asFieldOrder{k})        = v;
            end
        end
        fprintf("Frame-1 fallback: AS landmarks at AS index %d (abs frame %d).\n", ...
            asIdx_f1, AS_abs_start + asIdx_f1 - 1);
    end
end

addition = headingFromMVorPrev(seedMV1, validMask1, NaN, hasAS, states, startI, startI);

[params0, pHead0, torso0, seedAS] = buildParams0( ...
    seedMV1, validMask1, prev_pts_for_frame1, addition, ...
    scaleHead, scaleTorsoLength, Rx, Ry, Rz);

% Initialisation
f_seed = @(p) ASseedResidualFromModel( ...
    p, seedAS, Rz, Rx, Ry, scaleHead, scaleTorsoLength);

% Head AND torso angle windows are centred on their continuous seeds (head
% from the Kabsch fit, torso from the triangulated spine->neck axis), so the
% seed start is always inside the bounds and neither yaw is quantised to a
% 0/180 half-circle. theta_x_2 (torso roll, ~0) keeps the absolute +/-pi/5.
LB_seed = params0; UB_seed = params0;
LB_seed(1) = params0(1) - pi/5; UB_seed(1) = params0(1) + pi/5;
LB_seed(2) = params0(2) - pi/5; UB_seed(2) = params0(2) + pi/5;
LB_seed(3) = params0(3) - pi/2; UB_seed(3) = params0(3) + pi/2;
LB_seed(4:6) = pHead0.' - 0.5; UB_seed(4:6) = pHead0.' + 0.5;
LB_seed(7) = -pi/5;             UB_seed(7) = pi/5;
LB_seed(8) = params0(8) - pi/5; UB_seed(8) = params0(8) + pi/5;
LB_seed(9) = params0(9) - pi/2; UB_seed(9) = params0(9) + pi/2;
LB_seed(10:12) = torso0.' - 0.5; UB_seed(10:12) = torso0.' + 0.5;

seedOptions = optimoptions(@lsqnonlin, ...
    'Display','off', ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-10, ...
    'MaxFunctionEvaluations',5e4, ...
    'MaxIterations',5e3);

[params0, ~, residual_seed] = lsqnonlin(f_seed, params0, LB_seed, UB_seed, seedOptions);

fprintf('\n=== Frame 1 multiview-seed fit ===\n');
fprintf('Head init theta = [%.3f %.3f %.3f] deg\n', ...
    rad2deg(params0(1)), rad2deg(params0(2)), rad2deg(wrapToPi(params0(3))));
fprintf('Body init theta = [%.3f %.3f %.3f] deg\n', ...
    rad2deg(params0(7)), rad2deg(params0(8)), rad2deg(wrapToPi(params0(9))));
fprintf('Seed 3D fit residual norm = %.6g\n', norm(residual_seed));


% -------------------- LB/UB for MAIN optimizer (yaw widened to +/-pi) --------------------
xcams = zeros(1, ncams); ycams = zeros(1, ncams); zcams = zeros(1, ncams);
for ii = 1:ncams
    T = poses{ii}.Translation;
    xcams(ii) = T(1); ycams(ii) = T(2); zcams(ii) = T(3);
end

LB = zeros(1,12); UB = zeros(1,12);
LB(1)  = -pi/5;            UB(1)  =  pi/5;
LB(2)  = -pi/3;            UB(2)  =  pi/3;   % head pitch widened (it pinned at +/-pi/5 during the flip)
LB(3)  = -pi;              UB(3)  =  pi;
LB(4)  = min(xcams) - 5;   UB(4)  = max(xcams) + 5;
LB(5)  = min(ycams) - 2;   UB(5)  = max(ycams) + 2;
LB(6)  = min(zcams) - 0.5; UB(6)  = max(zcams) + 1.5;
LB(7)  = -pi/5;            UB(7)  =  pi/5;
LB(8)  = -pi/3;            UB(8)  =  pi/3;   % torso pitch widened (was +/-pi/5; avoid clipping a steeper body pitch)
LB(9)  = -pi;              UB(9)  =  pi;
LB(10) = min(xcams) - 5;   UB(10) = max(xcams) + 5;
LB(11) = min(ycams) - 2;   UB(11) = max(ycams) + 2;
LB(12) = min(zcams) - 0.5; UB(12) = max(zcams) + 1.5;

bounds = table(LB(:), UB(:), ...
    'VariableNames', {'Lower','Upper'}, ...
    'RowNames', ["theta_x","theta_y","theta_z","x","y","z", ...
                 "theta_x_2","theta_y_2","theta_z_2","x2","y2","z2"]);
disp(bounds);

% -------------------- LLH tensor --------------------
numvars = length(varsUsed);
LLH = zeros(nframes, numvars, ncams);
for ii = 1:ncams
    cam = camIDs(ii);
    camname = sprintf('cam%d', cam);
    camLLHmat = zeros(nframes, numvars);
    for j = 1:numvars
        bodypart = varsUsed{j};
        camLLHmat(:, j) = likelihoods.(camname).(bodypart);
    end
    LLH(:,:,ii) = camLLHmat;
end

% -------------------- Variances (ears included for CRBR) --------------------
clear variances;
variances.l_eye     = 3;
variances.r_eye     = 3;
variances.nose      = 2;
variances.l_ear     = 5;
variances.r_ear     = 5;
variances.spine     = 6;
variances.neck_base = 6;

% -------------------- Robust outlier handling in the fit --------------------
% Per-observation down-weighting so one inconsistent camera (e.g. a desynced /
% miscalibrated view sitting ~20px off) cannot drag the pose. Passed ONLY to
% the lsqnonlin solves below, NOT to the diagnostic reprojection-RMSE passes
% (those stay raw so the reported error is the true reprojection error).
%   type : 'none'  - plain least squares (cheetahError's fallback when robustOpt
%                    is omitted, used only by the diagnostic RMSE passes)
%          'huber' - residuals beyond c sigma get linear (not squared) influence
%          'tukey' - redescending: residuals beyond c sigma are smoothly ignored
%          'hard'  - residuals beyond c sigma are set to zero (rejected)
%   c    : threshold in SIGMA units. Residuals are (proj-label)/sigma, and good
%          observations already sit around ~2-3 sigma (the sigmas above are
%          optimistic), so use c~3 for 'huber' and c~4-6 for 'tukey'/'hard' to
%          avoid clipping good points while still removing a ~20px (>=8 sigma) view.
% Production reconstruction uses Huber (c=3) and passes it to every solve below;
% override via opts.robustOpt (e.g. struct('type','none') for OLS).
robustOpt = getOpt(opts, 'robustOpt', struct('type','huber','c',3));

% -------------------- Optimise frame 1 (useTemporal = false) --------------------
undistLabels = undistortedLabels;

f = @(params) cheetahError( ...
    params, Rz, Rx, Ry, scaleHead, scaleTorsoLength, ...
    undistLabels, poses, intrinsicParams, LLH, ...
    1, variances, params0, camIDs, false, robustOpt, camMaskStore(1,:));

options = optimoptions(@lsqnonlin, ...
    'Display','off', ...
    'StepTolerance',1e-4, ...
    'OptimalityTolerance',1e-4, ...
    'MaxFunctionEvaluations',1e5, ...
    'MaxIterations',1e5);

% Initial guess for the frame-1 solve (the multiview seed after the freeHead
% seed fit). Captured here so it can be compared against the optimised solution.
params0_frame1 = params0;

[xBar, ~, residual, ~, ~] = lsqnonlin(f, params0, LB, UB, options);

fprintf('\nOptimized frame 1 angles:\n');
fprintf('Head theta = [%.3f %.3f %.3f] deg\n', ...
    rad2deg(xBar(1)), rad2deg(xBar(2)), rad2deg(wrapToPi(xBar(3))));
fprintf('Body theta = [%.3f %.3f %.3f] deg\n', ...
    rad2deg(xBar(7)), rad2deg(xBar(8)), rad2deg(wrapToPi(xBar(9))));

% -------------------- Frame 1: initial guess vs solution vs AS (3D overlay) --------------------
% Single 3D axes overlaying the multiview-seeded INITIAL GUESS (params0_frame1,
% where the frame-1 solve started), the optimised SOLUTION (xBar), and the
% AcinoSet (AS) 3D landmarks for the same frame, so their agreement / offset can
% be inspected directly. All three live in the same world frame as the cameras.
% AS carries no ears, so its head skeleton uses eyes + nose only.
if show_frame1_init_vs_sol
    fprintf('\n=== Frame 1: initial guess vs solution ===\n');
    fprintf('                head theta [deg]              body theta [deg]\n');
    fprintf('init guess : [%7.2f %7.2f %7.2f]    [%7.2f %7.2f %7.2f]\n', ...
        rad2deg(params0_frame1(1)), rad2deg(params0_frame1(2)), rad2deg(wrapToPi(params0_frame1(3))), ...
        rad2deg(params0_frame1(7)), rad2deg(params0_frame1(8)), rad2deg(wrapToPi(params0_frame1(9))));
    fprintf('solution   : [%7.2f %7.2f %7.2f]    [%7.2f %7.2f %7.2f]\n', ...
        rad2deg(xBar(1)), rad2deg(xBar(2)), rad2deg(wrapToPi(xBar(3))), ...
        rad2deg(xBar(7)), rad2deg(xBar(8)), rad2deg(wrapToPi(xBar(9))));
    fprintf('init head pos = [%.3f %.3f %.3f], solution head pos = [%.3f %.3f %.3f]\n', ...
        params0_frame1(4), params0_frame1(5), params0_frame1(6), xBar(4), xBar(5), xBar(6));

    [~, ~, ptsInit] = cheetah3D(params0_frame1, Rz, Rx, Ry, scaleHead, scaleTorsoLength);
    [~, ~, ptsSol]  = cheetah3D(xBar,           Rz, Rx, Ry, scaleHead, scaleTorsoLength);

    lmNames   = {'l_eye','r_eye','nose','l_ear','r_ear','neck_base','spine'};
    lmColors  = {'b','r','k','c','m','g','y'};
    initColor = [0.6 0.6 0.6];
    asColor   = [0.95 0.45 0];

    figure('Color','w', ...
        'Name', sprintf('Frame 1 init / solution / AS 3D overlay (abs frame %d)', startI), ...
        'NumberTitle','off');
    set(gcf, 'Position', [120 80 1100 850]);
    hold on;

    % --- initial guess: gray hollow markers + dashed skeleton ---
    for k = 1:numel(lmNames)
        p = ptsInit.(lmNames{k})(:).';
        plot3(p(1), p(2), p(3), 'o', 'Color', initColor, 'MarkerSize', 6, 'LineWidth', 1.0);
    end
    drawSkeleton3D(ptsInit, initColor, '--', 1.5);

    % --- solution: per-landmark filled markers + solid black skeleton ---
    for k = 1:numel(lmNames)
        p = ptsSol.(lmNames{k})(:).';
        scatter3(p(1), p(2), p(3), 60, lmColors{k}, 'filled');
    end
    drawSkeleton3D(ptsSol, 'k', '-', 1.8);

    % --- AS (AcinoSet) seed: FULL landmark cloud (orange) + head/neck skeleton ---
    % Every AcinoSet landmark at this frame is plotted straight from AS.positions
    % (matching the complete source landmark cloud), not just the head/neck
    % subset, so the complete AcinoSet body is visible behind the model poses.
    % NaN rows (untracked landmarks) are dropped by scatter3 automatically.
    haveAS_f1 = hasAS && ~isnan(asSeedIdx);
    if haveAS_f1
        scatter3(AS.positions(:,1,asSeedIdx), AS.positions(:,2,asSeedIdx), ...
            AS.positions(:,3,asSeedIdx), 24, asColor, 'filled', 'o');
        % head/neck AS skeleton (eyes/nose, spine->neck) for orientation context
        drawSkeleton3D(asPts_frame1, asColor, ':', 1.6);
    else
        fprintf('Frame 1 overlay: no AS landmarks available (AS not loaded or FTE_start_frame.xlsx missing).\n');
    end

    % Legend proxies: one handle per pose (landmark colours apply to the solution).
    hInit = plot3(nan, nan, nan, 'o--', 'Color', initColor, 'LineWidth', 1.5, 'MarkerSize', 6);
    hSol  = plot3(nan, nan, nan, 'k-',  'LineWidth', 1.8);
    if haveAS_f1
        hAS = plot3(nan, nan, nan, 'o:', 'Color', asColor, ...
            'MarkerFaceColor', asColor, 'LineWidth', 1.6, 'MarkerSize', 7);
        legend([hInit hSol hAS], {'initial guess','solution','AS seed (full)'}, 'Location','best');
    else
        legend([hInit hSol], {'initial guess','solution'}, 'Location','best');
    end

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf(['Frame 1 (abs %d) 3D overlay: gray dashed = init guess, ' ...
        'colored solid = solution, orange = full AS seed'], startI));
    grid on; axis equal; view(3); rotate3d on; hold off;
end

% -------------------- Diagnose frame 1 --------------------
% CRBR frame-1 diagnostics: residual summary plus 3D scene,
% top/side/front ortho views, and a raw-video reprojection overlay of the
% optimised frame-1 solution. The diagnostic RMSE pass uses raw residuals
% (robustOpt is NOT applied) so the reported error is the true reprojection error.
if diagnose_frame_1
    % ----- Frame-1 residual quality check (sigma-normalized units) -----
    [E_f1, ~] = cheetahError( ...
        xBar, Rz, Rx, Ry, scaleHead, scaleTorsoLength, ...
        undistLabels, poses, intrinsicParams, LLH, ...
        1, variances, params0, camIDs, false);

    reprojRMSE_cam_f1 = zeros(1, ncams);
    R_all_f1 = [];
    for ii = 1:ncams
        R_cam = E_f1(:,:,ii);
        reprojRMSE_cam_f1(ii) = sqrt(mean(R_cam(:).^2));
        R_all_f1 = [R_all_f1; R_cam];
    end
    coverage_pairs_f1 = sum(any(R_all_f1 ~= 0, 2));
    reprojRMSE_f1     = sqrt(mean(R_all_f1(:).^2));

    fprintf("\n=== Frame 1 residual summary (sigma-normalized units) ===\n");
    fprintf("Frame 1 reprojection RMSE: %.6f\n", reprojRMSE_f1);
    fprintf("Frame 1 coverage pairs: %d\n", coverage_pairs_f1);
    for ii = 1:ncams
        fprintf("Frame 1 cam%d reprojection RMSE: %.6f\n", camIDs(ii), reprojRMSE_cam_f1(ii));
    end

    % ----- Reconstructed frame-1 keypoints and axes -----
    [~, ~, cheetahpts] = cheetah3D(xBar, Rz, Rx, Ry, scaleHead, scaleTorsoLength);

    l_eye     = cheetahpts.l_eye(:).';
    r_eye     = cheetahpts.r_eye(:).';
    nose      = cheetahpts.nose(:).';
    l_ear     = cheetahpts.l_ear(:).';
    r_ear     = cheetahpts.r_ear(:).';
    spine     = cheetahpts.spine(:).';
    neck_base = cheetahpts.neck_base(:).';

    coe      = (l_eye + r_eye) / 2;
    headAxis = nose - coe;         % centre-of-eyes -> nose
    bodyAxis = neck_base - spine;  % spine -> neck_base

    fprintf("\nReconstructed frame 1 (absolute frame %d)\n", startI);
    fprintf("Head axis (coe -> nose) = [%.6f %.6f %.6f]\n", ...
        headAxis(1), headAxis(2), headAxis(3));
    fprintf("Body axis (spine -> neck) = [%.6f %.6f %.6f]\n", ...
        bodyAxis(1), bodyAxis(2), bodyAxis(3));

    % Normalized arrows for plotting
    Lh = 0.30;
    Lb = 0.45;
    headAxisPlot = headAxis / max(norm(headAxis), 1e-12) * Lh;
    bodyAxisPlot = bodyAxis / max(norm(bodyAxis), 1e-12) * Lb;

    % ----- 3D scene plot -----
    figure('Color','w');
    set(gcf, 'Position', [80 80 1600 900]);
    if ncams == 6
        plotCameraConfig6(poses, 0.1);
    else
        plotCameraConfigLT6(poses, 0.1, camIDs);
    end
    hold on;

    scatter3(l_eye(1),     l_eye(2),     l_eye(3),     60, 'b', 'filled');
    scatter3(r_eye(1),     r_eye(2),     r_eye(3),     60, 'r', 'filled');
    scatter3(nose(1),      nose(2),      nose(3),      60, 'k', 'filled');
    scatter3(l_ear(1),     l_ear(2),     l_ear(3),     60, 'c', 'filled');
    scatter3(r_ear(1),     r_ear(2),     r_ear(3),     60, 'm', 'filled');
    scatter3(neck_base(1), neck_base(2), neck_base(3), 60, 'g', 'filled');
    scatter3(spine(1),     spine(2),     spine(3),     60, 'y', 'filled');

    plot3([l_eye(1) r_eye(1)], [l_eye(2) r_eye(2)], [l_eye(3) r_eye(3)], 'k-', 'LineWidth', 2);
    plot3([coe(1) nose(1)], [coe(2) nose(2)], [coe(3) nose(3)], 'k-', 'LineWidth', 2);
    plot3([neck_base(1) spine(1)], [neck_base(2) spine(2)], [neck_base(3) spine(3)], 'k-', 'LineWidth', 2);

    quiver3(coe(1), coe(2), coe(3), ...
            headAxisPlot(1), headAxisPlot(2), headAxisPlot(3), ...
            0, 'Color', [1 0 1], 'LineWidth', 3, 'MaxHeadSize', 1.4);
    quiver3(spine(1), spine(2), spine(3), ...
            bodyAxisPlot(1), bodyAxisPlot(2), bodyAxisPlot(3), ...
            0, 'Color', [0 0.2 1], 'LineWidth', 3, 'MaxHeadSize', 1.4);

    text(coe(1), coe(2), coe(3), '  coe', 'FontSize', 12, 'Color', [1 0 1]);
    text(nose(1), nose(2), nose(3), '  nose', 'FontSize', 12, 'Color', 'k');
    text(spine(1), spine(2), spine(3), '  spine', 'FontSize', 12, 'Color', [0 0.2 1]);
    text(neck_base(1), neck_base(2), neck_base(3), '  neck', 'FontSize', 12, 'Color', 'g');

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('Reconstructed 3D scene, frame 1 (absolute frame %d)', startI));
    grid on; axis equal; view(3); hold off;

    % ----- Top / Side / Front diagnostic views -----
    figure('Color','w');
    set(gcf, 'Position', [80 80 1800 600]);

    % Top view (X-Y)
    subplot(1,3,1); hold on;
    scatter(l_eye(1),     l_eye(2),     60, 'b', 'filled');
    scatter(r_eye(1),     r_eye(2),     60, 'r', 'filled');
    scatter(nose(1),      nose(2),      60, 'k', 'filled');
    scatter(l_ear(1),     l_ear(2),     60, 'c', 'filled');
    scatter(r_ear(1),     r_ear(2),     60, 'm', 'filled');
    scatter(neck_base(1), neck_base(2), 60, 'g', 'filled');
    scatter(spine(1),     spine(2),     60, 'y', 'filled');
    quiver(coe(1), coe(2), headAxisPlot(1), headAxisPlot(2), 0, ...
        'Color', [1 0 1], 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
    quiver(spine(1), spine(2), bodyAxisPlot(1), bodyAxisPlot(2), 0, ...
        'Color', [0 0.2 1], 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
    plot([l_eye(1) r_eye(1)], [l_eye(2) r_eye(2)], 'k-', 'LineWidth', 2);
    plot([coe(1) nose(1)], [coe(2) nose(2)], 'k-', 'LineWidth', 2);
    plot([neck_base(1) spine(1)], [neck_base(2) spine(2)], 'k-', 'LineWidth', 2);
    xlabel('X'); ylabel('Y'); title('Top view (X-Y)');
    axis equal; grid on;

    % Side view (X-Z)
    subplot(1,3,2); hold on;
    scatter(l_eye(1),     l_eye(3),     60, 'b', 'filled');
    scatter(r_eye(1),     r_eye(3),     60, 'r', 'filled');
    scatter(nose(1),      nose(3),      60, 'k', 'filled');
    scatter(l_ear(1),     l_ear(3),     60, 'c', 'filled');
    scatter(r_ear(1),     r_ear(3),     60, 'm', 'filled');
    scatter(neck_base(1), neck_base(3), 60, 'g', 'filled');
    scatter(spine(1),     spine(3),     60, 'y', 'filled');
    quiver(coe(1), coe(3), headAxisPlot(1), headAxisPlot(3), 0, ...
        'Color', [1 0 1], 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
    quiver(spine(1), spine(3), bodyAxisPlot(1), bodyAxisPlot(3), 0, ...
        'Color', [0 0.2 1], 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
    plot([l_eye(1) r_eye(1)], [l_eye(3) r_eye(3)], 'k-', 'LineWidth', 2);
    plot([coe(1) nose(1)], [coe(3) nose(3)], 'k-', 'LineWidth', 2);
    plot([neck_base(1) spine(1)], [neck_base(3) spine(3)], 'k-', 'LineWidth', 2);
    xlabel('X'); ylabel('Z'); title('Side view (X-Z)');
    axis equal; grid on;

    % Front view (Y-Z)
    subplot(1,3,3); hold on;
    scatter(l_eye(2),     l_eye(3),     60, 'b', 'filled');
    scatter(r_eye(2),     r_eye(3),     60, 'r', 'filled');
    scatter(nose(2),      nose(3),      60, 'k', 'filled');
    scatter(l_ear(2),     l_ear(3),     60, 'c', 'filled');
    scatter(r_ear(2),     r_ear(3),     60, 'm', 'filled');
    scatter(neck_base(2), neck_base(3), 60, 'g', 'filled');
    scatter(spine(2),     spine(3),     60, 'y', 'filled');
    quiver(coe(2), coe(3), headAxisPlot(2), headAxisPlot(3), 0, ...
        'Color', [1 0 1], 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
    quiver(spine(2), spine(3), bodyAxisPlot(2), bodyAxisPlot(3), 0, ...
        'Color', [0 0.2 1], 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
    plot([l_eye(2) r_eye(2)], [l_eye(3) r_eye(3)], 'k-', 'LineWidth', 2);
    plot([coe(2) nose(2)], [coe(3) nose(3)], 'k-', 'LineWidth', 2);
    plot([neck_base(2) spine(2)], [neck_base(3) spine(3)], 'k-', 'LineWidth', 2);
    xlabel('Y'); ylabel('Z'); title('Front view (Y-Z)');
    axis equal; grid on;

    % ----- Reproject reconstructed frame 1 to raw video -----
    cam = camIDs(1);   % choose a camera to overlay

    videoPath = fullfile(parent_folder, sprintf('cam%d_mjpeg.avi', cam));
    if exist(videoPath, 'file') ~= 2
        videoPath = fullfile(parent_folder, sprintf('cam%d.mp4', cam));
    end
    assert(exist(videoPath, 'file') == 2, 'Video not found: %s', videoPath);

    labelsToShow = {'l_eye','r_eye','nose','l_ear','r_ear','neck_base','spine'};

    camLocal = find(camIDs == cam, 1);
    assert(~isempty(camLocal), "Requested camera %d is not in camIDs.", cam);
    row    = find(offsets.cam == cam, 1);
    assert(~isempty(row), "Camera %d not found in offsets.", cam);
    offset = offsets.offset(row);

    frameIdx    = 1;   % xBar is the optimised result for frame 1 only
    rawFrameIdx = (startI + offset) + (frameIdx - 1);

    fprintf("Frame 1 reprojection: cam%d raw frame = %d (startI=%d, offset=%d)\n", ...
        cam, rawFrameIdx, startI, offset);

    vr = VideoReader(videoPath);
    t  = (rawFrameIdx - 1) / vr.FrameRate;
    t  = max(0, min(t, vr.Duration - 1/vr.FrameRate));
    vr.CurrentTime = t;
    img = readFrame(vr);

    % Build 3D points from optimised frame-1 solution
    [~, ~, pts] = cheetah3D(xBar, Rz, Rx, Ry, scaleHead, scaleTorsoLength);

    W = nan(numel(labelsToShow), 3);
    for k = 1:numel(labelsToShow)
        nm = labelsToShow{k};
        if isfield(pts, nm)
            p = pts.(nm); p = p(:)';
            if numel(p) == 3 && all(isfinite(p))
                W(k,:) = p;
            end
        end
    end

    validRecon = all(isfinite(W), 2);
    Wv         = W(validRecon, :);
    namesRecon = labelsToShow(validRecon);

    % Project world -> distorted image coordinates (matches raw video frame)
    uv = world2img(Wv, invert(poses{camLocal}), intrinsicParams{camLocal}, 'ApplyDistortion', true);

    % Load ORIGINAL RAW manual labels for this camera/frame
    camLabelsRaw = eval(sprintf('cam%dlabels', cam));
    assert(rawFrameIdx >= 1 && rawFrameIdx <= height(camLabelsRaw), ...
        'rawFrameIdx=%d outside cam%dlabels range 1..%d', ...
        rawFrameIdx, cam, height(camLabelsRaw));
    rawRow = camLabelsRaw(rawFrameIdx, :);

    uv_manual   = nan(numel(labelsToShow), 2);
    validManual = false(numel(labelsToShow), 1);
    for k = 1:numel(labelsToShow)
        bp    = labelsToShow{k};
        xName = sprintf('%s_x', bp);
        yName = sprintf('%s_y', bp);
        if ismember(xName, rawRow.Properties.VariableNames) && ...
           ismember(yName, rawRow.Properties.VariableNames)
            x = rawRow.(xName);
            y = rawRow.(yName);
            if ~isempty(x) && ~isempty(y) && all(isfinite([x y]))
                uv_manual(k,:) = [x y];
                validManual(k) = true;
            end
        end
    end

    % Overlay raw labels + reprojection
    figure('Color','w', ...
        'Name', sprintf('Cam%d reprojection: frame 1 -> raw frame %d', cam, rawFrameIdx), ...
        'NumberTitle', 'off');
    imshow(img); hold on;

    plot(uv(:,1), uv(:,2), 'go', 'MarkerSize', 8, 'LineWidth', 2);
    for k = 1:numel(namesRecon)
        text(uv(k,1), uv(k,2), ['  reproj-' namesRecon{k}], ...
            'Color','g', 'FontSize', 11, 'FontWeight','bold');
    end

    plot(uv_manual(validManual,1), uv_manual(validManual,2), ...
        'r+', 'MarkerSize', 10, 'LineWidth', 2);
    for k = 1:numel(labelsToShow)
        if validManual(k)
            text(uv_manual(k,1), uv_manual(k,2), ['  manual-' labelsToShow{k}], ...
                'Color','r', 'FontSize', 11, 'FontWeight','bold');
        end
    end

    legend({'reprojected 3D','manual raw labels'}, 'TextColor','w', 'Location','best');
    title(sprintf('Cam%d raw frame %d: green = reprojection, red = manual labels', ...
        cam, rawFrameIdx));
    hold off;
end

% -------------------- Main optimisation loop (frames 2:nframes) --------------------
% Per-frame multiview-seeded init, useTemporal = false on every frame.

residualStore     = residual;
xBarStore         = xBar;
seedSourceLog     = strings(nframes, 1);
nValidMVperFrame  = zeros(nframes, 1);

seedSourceLog(1)    = "MV_frame1";
nValidMVperFrame(1) = nValid1;

seedRefineOpts = optimoptions(@lsqnonlin, ...
    'Display','off', ...
    'StepTolerance',1e-8, ...
    'OptimalityTolerance',1e-8, ...
    'MaxFunctionEvaluations',1e4, ...
    'MaxIterations',1e3);

for i = 2:nframes
    abs_frame_i = startI + (i - 1);
    prevSol     = xBarStore(i-1, :);
    [~, ~, prev_pts] = cheetah3D(prevSol, Rz, Rx, Ry, scaleHead, scaleTorsoLength);

    [seedMV_i, validMask_i] = getMVseed(abs_frame_i, triRowOf, TRI, landmarkNames);
    nValid_i = sum(validMask_i);
    nValidMVperFrame(i) = nValid_i;

    if nValid_i < 2
        % Too sparse to seed independently; reuse previous solution.
        initGuess        = prevSol;
        seedSourceLog(i) = "prev_pose_sparse_MV";
    else
        addition_i = headingFromMVorPrev( ...
            seedMV_i, validMask_i, prevSol(9), hasAS, states, abs_frame_i, startI);

        [params0_i, pHead0_i, torso0_i, sf_i] = buildParams0( ...
            seedMV_i, validMask_i, prev_pts, addition_i, ...
            scaleHead, scaleTorsoLength, Rx, Ry, Rz);

        % Head AND torso windows centred on their continuous seeds (head Kabsch,
        % torso triangulated spine->neck axis); no 0/180 yaw quantisation.
        % theta_x_2 (torso roll, ~0) keeps the absolute +/-pi/5.
        LB_seed_i = params0_i; UB_seed_i = params0_i;
        LB_seed_i(1) = params0_i(1) - pi/5; UB_seed_i(1) = params0_i(1) + pi/5;
        LB_seed_i(2) = params0_i(2) - pi/5; UB_seed_i(2) = params0_i(2) + pi/5;
        LB_seed_i(3) = params0_i(3) - pi/2; UB_seed_i(3) = params0_i(3) + pi/2;
        LB_seed_i(4:6) = pHead0_i.' - 0.5; UB_seed_i(4:6) = pHead0_i.' + 0.5;
        LB_seed_i(7) = -pi/5;               UB_seed_i(7) = pi/5;
        LB_seed_i(8) = params0_i(8) - pi/5; UB_seed_i(8) = params0_i(8) + pi/5;
        LB_seed_i(9) = params0_i(9) - pi/2; UB_seed_i(9) = params0_i(9) + pi/2;
        LB_seed_i(10:12) = torso0_i.' - 0.5; UB_seed_i(10:12) = torso0_i.' + 0.5;

        f_seed_i = @(p) ASseedResidualFromModel( ...
            p, sf_i, Rz, Rx, Ry, scaleHead, scaleTorsoLength);

        try
            params0_i = lsqnonlin(f_seed_i, params0_i, LB_seed_i, UB_seed_i, seedRefineOpts);
            if nValid_i >= 4
                seedSourceLog(i) = "MV_refined";
            else
                seedSourceLog(i) = "MV_partial";
            end
        catch
            seedSourceLog(i) = "prev_pose_seed_fail";
            params0_i = prevSol;
        end

        initGuess = params0_i;
    end

    % Clamp into main bounds (per-frame seed may sit at bound edges)
    initGuess = min(max(initGuess, LB), UB);

    % Main lsqnonlin: reprojection-only, no temporal coupling. camMaskStore(i,:)
    % hard-excludes any camera the LOO trust pass flagged as biased this frame.
    f = @(params) cheetahError( ...
        params, Rz, Rx, Ry, scaleHead, scaleTorsoLength, ...
        undistLabels, poses, intrinsicParams, LLH, ...
        i, variances, prevSol, camIDs, false, robustOpt, camMaskStore(i,:));

    % --- Multi-start keep-best: continuous Kabsch seed vs previous-frame pose ---
    % The per-frame solve can fall into a flipped-head local minimum when the
    % triangulated head seed is momentarily ambiguous (e.g. mid-turn). Solving
    % from both the fresh seed and the previous good pose, then keeping the lower
    % RAW-reprojection result, prevents the solver from getting trapped.
    prevStart = min(max(prevSol, LB), UB);
    starts    = {initGuess, prevStart};

    % Select on the ROBUST (consensus-rejected) cost the optimiser actually
    % minimised, NOT raw reprojection RMSE. A systematically biased camera (e.g. a
    % close, mis-calibrated / sub-frame-desynced view) makes a flipped-head pose
    % look good in RAW pixels, so raw selection can pick the wrong basin. The
    % robust residual excludes the rejected outlier views, so the pose consistent
    % with the camera majority wins.
    bestXBar = []; bestRes = []; bestCost = inf; bestStart = 1;
    for sIdx = 1:numel(starts)
        [xCand, resnormCand, resCand, ~, ~] = lsqnonlin(f, starts{sIdx}, LB, UB, options);
        if resnormCand < bestCost
            bestCost = resnormCand; bestXBar = xCand; bestRes = resCand; bestStart = sIdx;
        end
    end
    xBar     = bestXBar;
    residual = bestRes;
    if bestStart == 2
        seedSourceLog(i) = seedSourceLog(i) + "+prevWin";
    end

    residualStore = vertcat(residualStore, residual);
    xBarStore     = vertcat(xBarStore, xBar);
end

% -------------------- Seed source summary --------------------
fprintf("\n=== Per-frame seed source summary ===\n");
uniqueSources = unique(seedSourceLog);
for s = 1:numel(uniqueSources)
    src = uniqueSources(s);
    n   = sum(seedSourceLog == src);
    fprintf("  %-22s : %d frames\n", src, n);
end
fprintf("nValidMV per frame: mean=%.2f, median=%d, min=%d, max=%d\n", ...
    mean(nValidMVperFrame), median(nValidMVperFrame), min(nValidMVperFrame), max(nValidMVperFrame));

if doFigs
figure('Color','w','Name','Multiview landmark coverage per frame', 'Visible', onoff(showFigures));
plot(1:nframes, nValidMVperFrame, 'k-', 'LineWidth', 1.2);
yline(7, 'g--', 'full'); yline(4, 'b--', 'min ok'); yline(2, 'r--', 'fallback');
grid on; xlabel('Frame'); ylabel('# valid MV landmarks');
title('Per-frame multiview triangulation coverage');
xlim([1 nframes]);
if saveFigures
    print(gcf, fullfile(parent_folder, resultsSubdir, 'MV_coverage_per_frame.png'), '-dpng', '-r200');
end
end

% -------------------- Residual quality checks (per frame, reprojection only) --------------------
reprojRMSE     = zeros(nframes, 1);
coverage_pairs = zeros(nframes, 1);
reprojRMSE_cam = zeros(nframes, ncams);

for fidx = 1:nframes
    [E, ~] = cheetahError( ...
        xBarStore(fidx,:), Rz, Rx, Ry, scaleHead, scaleTorsoLength, ...
        undistLabels, poses, intrinsicParams, LLH, ...
        fidx, variances, xBarStore(fidx,:), camIDs, false);

    R_all = [];
    for cam = 1:ncams
        R_cam = E(:,:,cam);
        reprojRMSE_cam(fidx, cam) = sqrt(mean(R_cam(:).^2));
        R_all = [R_all; R_cam];
    end
    coverage_pairs(fidx) = sum(any(R_all ~= 0, 2));
    reprojRMSE(fidx)     = sqrt(mean(R_all(:).^2));
end

fprintf("\n=== Residual summary (sigma-normalized units, reprojection only) ===\n");
fprintf("Reproj RMSE: mean=%.4f, median=%.4f, p95=%.4f, max=%.4f\n", ...
    mean(reprojRMSE), median(reprojRMSE), prctile(reprojRMSE,95), max(reprojRMSE));
fprintf("Coverage pairs per frame: mean=%.2f, min=%d, max=%d\n", ...
    mean(coverage_pairs), min(coverage_pairs), max(coverage_pairs));

if doFigs
figure('Visible', onoff(showFigures)); plot(reprojRMSE); grid on;
title("Per-frame reprojection RMSE (sigma-normalized)"); xlabel("Frame"); ylabel("RMSE");
if saveFigures
    print(gcf, fullfile(parent_folder, resultsSubdir, 'Per-frame_reprojection_RMSE.png'), '-dpng', '-r200');
end

figure('Visible', onoff(showFigures)); plot(coverage_pairs); grid on;
title("Per-frame observation coverage"); xlabel("Frame"); ylabel("# pairs");
if saveFigures
    print(gcf, fullfile(parent_folder, resultsSubdir, 'Per-frame_observation_coverage.png'), '-dpng', '-r200');
end

figure('Visible', onoff(showFigures)); plot(reprojRMSE_cam); grid on;
title("Per-camera reprojection RMSE (sigma-normalized)"); xlabel("Frame"); ylabel("RMSE");
legend(arrayfun(@(c) sprintf("cam%d", c), camIDs, "UniformOutput", false), "Location", "best");
if saveFigures
    print(gcf, fullfile(parent_folder, resultsSubdir, 'Per-camera_reprojection_RMSE.png'), '-dpng', '-r200');
end
end

% -------------------- Camera-trust transparency log (fit stage) --------------------
% Per-frame record of which cameras the leave-one-out trust pass kept vs
% hard-excluded from the cheetahError fit, with each camera's LOO residual (from
% triangulation) and fit reprojection RMSE, so the rejection is fully auditable.
nKeptPerFrame = sum(camMaskStore, 2);
absFrameCol   = (startI + (0:nframes-1)).';
trustTbl = table((1:nframes).', absFrameCol, nKeptPerFrame, ...
    'VariableNames', {'frame','absFrame','nKept'});
for c = 1:ncams
    cid = camIDs(c);
    trustTbl.(sprintf('cam%d_kept',       cid)) = double(camMaskStore(:,c));
    trustTbl.(sprintf('cam%d_LOOresid',   cid)) = camLOOStore(:,c);
    trustTbl.(sprintf('cam%d_reprojRMSE', cid)) = reprojRMSE_cam(:,c);
end
trustCsv = fullfile(parent_folder, resultsSubdir, 'camRejectionLog.csv');
writetable(trustTbl, trustCsv);
fprintf("\n=== Camera-trust (fit stage) summary ===\n");
fprintf("Per-frame camera rejection log: %s\n", trustCsv);
for c = 1:ncams
    fprintf("  cam%d: excluded from fit in %d / %d frames\n", ...
        camIDs(c), sum(~camMaskStore(:,c)), nframes);
end

if doFigs
figure('Color','w', 'Visible', onoff(showFigures)); hold on;
mvCamColors = lines(ncams);
for c = 1:ncams
    plot(1:nframes, reprojRMSE_cam(:,c), '-', 'Color', mvCamColors(c,:), 'LineWidth', 1.2);
end
for c = 1:ncams
    rej = find(~camMaskStore(:,c));
    if ~isempty(rej)
        plot(rej, reprojRMSE_cam(rej,c), 'x', 'Color', mvCamColors(c,:), ...
            'MarkerSize', 9, 'LineWidth', 1.6);
    end
end
grid on; xlabel("Frame"); ylabel("Reproj RMSE (sigma-norm)");
title("Per-camera reprojection RMSE (x = excluded by leave-one-out trust)");
legend(arrayfun(@(c) sprintf("cam%d", c), camIDs, "UniformOutput", false), "Location", "best");
if saveFigures
    print(gcf, fullfile(parent_folder, resultsSubdir, 'Per-camera_reprojection_RMSE_rejections.png'), '-dpng', '-r200');
end
end

% -------------------- Common sizes --------------------
sz       = size(xBarStore, 1);
frameIdx = 1:sz;

% -------------------- Median filter --------------------
window = max(round(sz * 0.0625), 3);

xBarFiltStore_Med = [ ...
    medfilt1(xBarStore(:,1),  window,   'truncate') ...
    medfilt1(xBarStore(:,2),  window,   'truncate') ...
    medfilt1(xBarStore(:,3),  window,   'truncate') ...
    medfilt1(xBarStore(:,4),  window,   'truncate') ...
    medfilt1(xBarStore(:,5),  window,   'truncate') ...
    medfilt1(xBarStore(:,6),  window,   'truncate') ...
    medfilt1(xBarStore(:,7),  window,   'truncate') ...
    medfilt1(xBarStore(:,8),  window,   'truncate') ...
    medfilt1(xBarStore(:,9),  window,   'truncate') ...
    medfilt1(xBarStore(:,10), window+3, 'truncate') ...
    medfilt1(xBarStore(:,11), window+3, 'truncate') ...
    medfilt1(xBarStore(:,12), window+3, 'truncate')];

save(fullfile(parent_folder, "storedOptimMed" + optimSuffix + ".mat"), ...
    "xBarFiltStore_Med", "xBarStore", "absFrames");

% -------------------- Low-pass filter --------------------
% fs is the dataset frame rate (90 or 120 Hz); runFullPipeline passes the
% detected rate via opts.fs. Fall back to 90 Hz for standalone calls, matching
% the repo-wide default. This filter feeds only the optional LPF trajectory.
fc = getOpt(opts, 'fc', 5);
fs = getOpt(opts, 'fs', 90);
[b, a] = butter(4, fc / (fs/2), 'low');

xBarFiltStore_LPF = [ ...
    filtfilt(b, a, xBarStore(:,1)) ...
    filtfilt(b, a, xBarStore(:,2)) ...
    filtfilt(b, a, xBarStore(:,3)) ...
    filtfilt(b, a, xBarStore(:,4)) ...
    filtfilt(b, a, xBarStore(:,5)) ...
    filtfilt(b, a, xBarStore(:,6)) ...
    filtfilt(b, a, xBarStore(:,7)) ...
    filtfilt(b, a, xBarStore(:,8)) ...
    filtfilt(b, a, xBarStore(:,9)) ...
    filtfilt(b, a, xBarStore(:,10)) ...
    filtfilt(b, a, xBarStore(:,11)) ...
    filtfilt(b, a, xBarStore(:,12))];

save(fullfile(parent_folder, "storedOptimLPF" + optimSuffix + ".mat"), ...
    "xBarFiltStore_LPF", "xBarStore", "absFrames");

% -------------------- Raw results --------------------
xBarRawStore = xBarStore;
% Provenance so downstream consumers (Monte Carlo, pursuit) can assert the
% trajectory is the frame-wise raw optimiser output produced with Huber.
reconMeta = struct('reconMode', "mvperframe", 'source', "raw", ...
    'temporal', false, 'robustType', string(robustOpt.type), ...
    'robustC', robustOpt.c, 'fs', fs, 'fc', fc, ...
    'absFrameStart', absFrames(1), 'absFrameEnd', absFrames(end), ...
    'nFrames', numel(absFrames));
save(fullfile(parent_folder, "storedOptimRaw" + optimSuffix + ".mat"), ...
    "xBarRawStore", "xBarStore", "seedSourceLog", "nValidMVperFrame", ...
    "absFrames", "reconMeta");

% -------------------- Choose trajectory source for video --------------------
switch lower(videoSource)
    case "raw"; xBarVideo = xBarRawStore;       videoTag = "raw";
    case "med"; xBarVideo = xBarFiltStore_Med;  videoTag = "med";
    case "lpf"; xBarVideo = xBarFiltStore_LPF;  videoTag = "lpf";
    otherwise
        error('videoSource must be "raw", "med", or "lpf"');
end

% -------------------- Orthographic top/side/front video --------------------
if createOrthoVideo
    writeName = fullfile(parent_folder, resultsSubdir, "LSQrecon_ortho_" + videoTag);
    myWriter = VideoWriter(writeName, 'Motion JPEG AVI');
    myWriter.Quality   = 100;
    myWriter.FrameRate = 20;
    open(myWriter);

    fg = figure('Color','w','Visible',onoff(showVideos)); set(fg, 'Position', [80 80 1800 600]);

    for i = 1:size(xBarVideo, 1)
        clf(fg);

        [~, ~, cheetahpts] = cheetah3D(xBarVideo(i,:), Rz, Rx, Ry, scaleHead, scaleTorsoLength);

        l_eye = cheetahpts.l_eye(:).'; r_eye = cheetahpts.r_eye(:).';
        nose  = cheetahpts.nose(:).';  neck_base = cheetahpts.neck_base(:).';
        l_ear = cheetahpts.l_ear(:).'; r_ear = cheetahpts.r_ear(:).';
        spine = cheetahpts.spine(:).';

        coe      = (l_eye + r_eye) / 2;
        headAxis = nose - coe;
        bodyAxis = neck_base - spine;

        Lh = 0.30; Lb = 0.45;
        headAxisPlot = headAxis / max(norm(headAxis), 1e-12) * Lh;
        bodyAxisPlot = bodyAxis / max(norm(bodyAxis), 1e-12) * Lb;

        allPts = [l_eye; r_eye; nose; l_ear; r_ear; neck_base; spine];
        xs = allPts(:,1); ys = allPts(:,2); zs = allPts(:,3);
        padXY = 0.08; padZ = 0.08;
        xlims_i = [min(xs)-padXY, max(xs)+padXY];
        ylims_i = [min(ys)-padXY, max(ys)+padXY];
        zlims_i = [min(zs)-padZ , max(zs)+padZ ];

        subplot(1,3,1); hold on;
        scatter(l_eye(1), l_eye(2), 60,'b','filled');
        scatter(r_eye(1), r_eye(2), 60,'r','filled');
        scatter(nose(1),  nose(2),  60,'k','filled');
        scatter(l_ear(1), l_ear(2), 60,'c','filled');
        scatter(r_ear(1), r_ear(2), 60,'m','filled');
        scatter(neck_base(1), neck_base(2), 60,'g','filled');
        scatter(spine(1), spine(2), 60,'y','filled');
        plot([l_eye(1) r_eye(1)],[l_eye(2) r_eye(2)],'k-','LineWidth',2);
        plot([coe(1) nose(1)],[coe(2) nose(2)],'k-','LineWidth',2);
        plot([l_ear(1) r_ear(1)],[l_ear(2) r_ear(2)],'k-','LineWidth',2);
        plot([neck_base(1) spine(1)],[neck_base(2) spine(2)],'k-','LineWidth',2);
        quiver(coe(1), coe(2), headAxisPlot(1), headAxisPlot(2), 0, ...
            'Color',[1 0 1],'LineWidth',2.5,'MaxHeadSize',1.5);
        quiver(spine(1), spine(2), bodyAxisPlot(1), bodyAxisPlot(2), 0, ...
            'Color',[0 0.2 1],'LineWidth',2.5,'MaxHeadSize',1.5);
        xlabel('X'); ylabel('Y');
        title(sprintf('Top view (X-Y), frame %d', i));
        axis equal; grid on; xlim(xlims_i); ylim(ylims_i);

        subplot(1,3,2); hold on;
        scatter(l_eye(1), l_eye(3), 60,'b','filled');
        scatter(r_eye(1), r_eye(3), 60,'r','filled');
        scatter(nose(1),  nose(3),  60,'k','filled');
        scatter(l_ear(1), l_ear(3), 60,'c','filled');
        scatter(r_ear(1), r_ear(3), 60,'m','filled');
        scatter(neck_base(1), neck_base(3), 60,'g','filled');
        scatter(spine(1), spine(3), 60,'y','filled');
        plot([l_eye(1) r_eye(1)],[l_eye(3) r_eye(3)],'k-','LineWidth',2);
        plot([coe(1) nose(1)],[coe(3) nose(3)],'k-','LineWidth',2);
        plot([l_ear(1) r_ear(1)],[l_ear(3) r_ear(3)],'k-','LineWidth',2);
        plot([neck_base(1) spine(1)],[neck_base(3) spine(3)],'k-','LineWidth',2);
        quiver(coe(1), coe(3), headAxisPlot(1), headAxisPlot(3), 0, ...
            'Color',[1 0 1],'LineWidth',2.5,'MaxHeadSize',1.5);
        quiver(spine(1), spine(3), bodyAxisPlot(1), bodyAxisPlot(3), 0, ...
            'Color',[0 0.2 1],'LineWidth',2.5,'MaxHeadSize',1.5);
        xlabel('X'); ylabel('Z');
        title(sprintf('Side view (X-Z), frame %d', i));
        axis equal; grid on; xlim(xlims_i); ylim(zlims_i);

        subplot(1,3,3); hold on;
        scatter(l_eye(2), l_eye(3), 60,'b','filled');
        scatter(r_eye(2), r_eye(3), 60,'r','filled');
        scatter(nose(2),  nose(3),  60,'k','filled');
        scatter(l_ear(2), l_ear(3), 60,'c','filled');
        scatter(r_ear(2), r_ear(3), 60,'m','filled');
        scatter(neck_base(2), neck_base(3), 60,'g','filled');
        scatter(spine(2), spine(3), 60,'y','filled');
        plot([l_eye(2) r_eye(2)],[l_eye(3) r_eye(3)],'k-','LineWidth',2);
        plot([coe(2) nose(2)],[coe(3) nose(3)],'k-','LineWidth',2);
        plot([l_ear(2) r_ear(2)],[l_ear(3) r_ear(3)],'k-','LineWidth',2);
        plot([neck_base(2) spine(2)],[neck_base(3) spine(3)],'k-','LineWidth',2);
        quiver(coe(2), coe(3), headAxisPlot(2), headAxisPlot(3), 0, ...
            'Color',[1 0 1],'LineWidth',2.5,'MaxHeadSize',1.5);
        quiver(spine(2), spine(3), bodyAxisPlot(2), bodyAxisPlot(3), 0, ...
            'Color',[0 0.2 1],'LineWidth',2.5,'MaxHeadSize',1.5);
        xlabel('Y'); ylabel('Z');
        title(sprintf('Front view (Y-Z), frame %d', i));
        axis equal; grid on; xlim(ylims_i); ylim(zlims_i);

        sgtitle(sprintf('Ortho views (CRBR MV per-frame, %s), abs frame %d', videoTag, startI + i - 1));

        drawnow;
        frame = getframe(fg);
        writeVideo(myWriter, frame);
    end

    close(myWriter); close(fg);
    fprintf('Orthographic video saved: %s.avi\n', writeName);
end

% -------------------- 3D scene recon video --------------------
if createVideo
    writeName = fullfile(parent_folder, resultsSubdir, "LSQrecon_" + videoTag);
    myWriter = VideoWriter(writeName, 'Motion JPEG AVI');
    myWriter.Quality   = 100;
    myWriter.FrameRate = 20;
    open(myWriter);

    fg = figure('Visible', onoff(showVideos));
    set(fg, 'Position', [0 0 1920 1080]);
    set(fg, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [1920 1080]);

    for i = 1:sz
        clf(fg);
        if ncams == 6
            plotCameraConfig6(poses, 0.1);
        else
            plotCameraConfigLT6(poses, 0.1, camIDs);
        end
        hold on;
        [~, ~, cheetahpts] = cheetah3D(xBarVideo(i,:), Rz, Rx, Ry, scaleHead, scaleTorsoLength);
        plotCheetah(cheetahpts); hold on;
        xlabel('x'); ylabel('y'); zlabel('z'); hold off;
        sgtitle('Reconstructed 3D Scene (CRBR MV per-frame, useTemporal=false)');
        frame = getframe(fg);
        writeVideo(myWriter, frame);
    end
    close(myWriter);

    outFile = writeName;
    if isstring(outFile), outFile = char(outFile); end
    if ~endsWith(lower(outFile), ".avi")
        outFile = [outFile '.avi'];
    end
    fprintf("Wrote 3D scene video: %s\n", outFile);
end

% ==================== VALIDATION: Reprojection across all frames & cameras ====================
fprintf('\n=== Starting full reprojection validation ===\n');

varnames_val = {'l_eye','r_eye','nose','l_ear','r_ear','spine','neck_base'};
nparts       = numel(varnames_val);

partColors = [
    0.0  0.4  1.0;   % l_eye     - blue
    1.0  0.2  0.2;   % r_eye     - red
    0.0  0.0  0.0;   % nose      - black
    0.0  0.8  0.8;   % l_ear     - cyan
    0.8  0.0  0.8;   % r_ear     - magenta
    0.9  0.7  0.0;   % spine     - yellow
    0.0  0.7  0.0;   % neck_base - green
];

val_errPix           = nan(nframes, ncams, nparts);
val_reprojRMSE_frame = nan(nframes, 1);
val_reprojRMSE_cam   = nan(nframes, ncams);
val_reprojRMSE_part  = nan(nframes, nparts);

% Mask-aware "kept cameras only" counterparts. Same reprojection error, but per
% frame they exclude any camera the leave-one-out trust pass (camMaskStore)
% rejected from that frame's fit. The all-camera arrays above still report the
% TRUE reprojection error against every visible label (so a biased camera is
% visible in the diagnostics); these report the error of just the views the pose
% was actually fit to, so one systematically biased camera no longer inflates the
% headline RMSE. camMaskStore is (nframes x ncams), rows aligned to fidx.
val_errPix_kept           = nan(nframes, ncams, nparts);
val_reprojRMSE_frame_kept = nan(nframes, 1);
val_reprojRMSE_cam_kept   = nan(nframes, ncams);
val_reprojRMSE_part_kept  = nan(nframes, nparts);

for fidx = 1:nframes
    [~, ~, cheetahpts_val] = cheetah3D( ...
        xBarStore(fidx,:), Rz, Rx, Ry, scaleHead, scaleTorsoLength);

    errAll_frame      = [];
    errAll_frame_kept = [];   % all-cam errors minus cameras rejected from this frame's fit

    for ii = 1:ncams
        cam     = camIDs(ii);
        camname = sprintf('cam%d', cam);
        camKept = camMaskStore(fidx, ii);   % false => LOO trust pass excluded this cam this frame
        errCam     = [];
        errCamKept = [];

        for j = 1:nparts
            bodypart    = varnames_val{j};
            coord_3d    = cheetahpts_val.(bodypart)';
            pixels_proj = world2img(coord_3d, invert(poses{ii}), intrinsicParams{ii});
            pixels_obs  = undistortedLabels.(camname).(bodypart){fidx};

            llhIdx = find(strcmp(varsUsed, bodypart), 1);
            valid  = ~isempty(llhIdx) && ...
                     LLH(fidx, llhIdx, ii) > 0.5 && ...
                     all(isfinite(pixels_proj)) && all(isfinite(pixels_obs)) && ...
                     all(pixels_proj > 1) && all(pixels_obs > 1);

            if valid
                errEuc = norm(pixels_proj - pixels_obs);
                val_errPix(fidx, ii, j) = errEuc;
                errCam       = [errCam;       errEuc];
                errAll_frame = [errAll_frame; errEuc];
                if camKept
                    val_errPix_kept(fidx, ii, j) = errEuc;
                    errCamKept        = [errCamKept;        errEuc];
                    errAll_frame_kept = [errAll_frame_kept; errEuc];
                end
            end
        end

        if ~isempty(errCam)
            val_reprojRMSE_cam(fidx, ii) = sqrt(mean(errCam.^2));
        end
        if ~isempty(errCamKept)
            val_reprojRMSE_cam_kept(fidx, ii) = sqrt(mean(errCamKept.^2));
        end
    end

    for j = 1:nparts
        partVals = squeeze(val_errPix(fidx, :, j));
        partVals = partVals(~isnan(partVals));
        if ~isempty(partVals)
            val_reprojRMSE_part(fidx, j) = sqrt(mean(partVals.^2));
        end

        partValsKept = squeeze(val_errPix_kept(fidx, :, j));
        partValsKept = partValsKept(~isnan(partValsKept));
        if ~isempty(partValsKept)
            val_reprojRMSE_part_kept(fidx, j) = sqrt(mean(partValsKept.^2));
        end
    end

    if ~isempty(errAll_frame)
        val_reprojRMSE_frame(fidx) = sqrt(mean(errAll_frame.^2));
    end
    if ~isempty(errAll_frame_kept)
        val_reprojRMSE_frame_kept(fidx) = sqrt(mean(errAll_frame_kept.^2));
    end
end

save(fullfile(parent_folder, "cheetahReprojErr" + optimSuffix + ".mat"), ...
    "val_reprojRMSE_frame", "val_reprojRMSE_cam", "val_reprojRMSE_part", ...
    "val_reprojRMSE_frame_kept", "val_reprojRMSE_cam_kept", "val_reprojRMSE_part_kept", ...
    "val_errPix", "val_errPix_kept", "varnames_val", "seedSourceLog", ...
    "nValidMVperFrame", "absFrames");

fprintf('\n--- Overall ---\n');
fprintf('Per-frame RMSE (all cams) : mean=%.3f median=%.3f p95=%.3f max=%.3f (px)\n', ...
    mean(val_reprojRMSE_frame,'omitnan'), ...
    median(val_reprojRMSE_frame,'omitnan'), ...
    prctile(val_reprojRMSE_frame(~isnan(val_reprojRMSE_frame)), 95), ...
    max(val_reprojRMSE_frame,[],'omitnan'));
vfk = val_reprojRMSE_frame_kept(~isnan(val_reprojRMSE_frame_kept));
fprintf('Per-frame RMSE (kept cams): mean=%.3f median=%.3f p95=%.3f max=%.3f (px)\n', ...
    mean(val_reprojRMSE_frame_kept,'omitnan'), ...
    median(val_reprojRMSE_frame_kept,'omitnan'), ...
    prctile(vfk, 95), ...
    max(val_reprojRMSE_frame_kept,[],'omitnan'));

fprintf('\n--- Per camera (all frames vs. kept-only) ---\n');
for ii = 1:ncams
    col  = val_reprojRMSE_cam(:,ii);
    colK = val_reprojRMSE_cam_kept(:,ii);
    nKept = sum(~isnan(colK));
    fprintf('cam%d: all  mean=%.3f median=%.3f p95=%.3f max=%.3f\n', ...
        camIDs(ii), mean(col,'omitnan'), median(col,'omitnan'), ...
        prctile(col(~isnan(col)),95), max(col,[],'omitnan'));
    if nKept > 0
        fprintf('      kept mean=%.3f median=%.3f p95=%.3f max=%.3f  (%d/%d frames kept)\n', ...
            mean(colK,'omitnan'), median(colK,'omitnan'), ...
            prctile(colK(~isnan(colK)),95), max(colK,[],'omitnan'), nKept, nframes);
    else
        fprintf('      kept: rejected in all frames with valid observations\n');
    end
end

fprintf('\n--- Per body part ---\n');
for j = 1:nparts
    allErr = squeeze(val_errPix(:,:,j));
    allErr = allErr(~isnan(allErr));
    if ~isempty(allErr)
        fprintf('%-12s mean=%.3f median=%.3f p95=%.3f max=%.3f\n', ...
            varnames_val{j}, mean(allErr), median(allErr), prctile(allErr,95), max(allErr));
    else
        fprintf('%-12s no valid observations\n', varnames_val{j});
    end
end

% -------------------- Validation figures --------------------
if doFigs
figure('Color','w','Name','Validation: Per-frame RMSE overview', 'Visible', onoff(showFigures), 'Position',[100 100 1200 450]);

subplot(1,2,1);
hold on;
plot(1:nframes, val_reprojRMSE_frame,      'k-', 'LineWidth', 1.5);
plot(1:nframes, val_reprojRMSE_frame_kept, '-',  'Color',[0 0.6 0], 'LineWidth', 1.5);
yline(mean(val_reprojRMSE_frame,'omitnan'),      'k--', 'LineWidth', 1.0);
yline(mean(val_reprojRMSE_frame_kept,'omitnan'), '--', 'Color',[0 0.6 0], 'LineWidth', 1.0);
grid on; xlabel('Frame'); ylabel('RMSE (px)');
title('Overall per-frame reprojection RMSE');
legend({'all cams','kept cams (LOO-masked)','all-cam mean','kept-cam mean'}, 'Location','best');
xlim([1 nframes]);

subplot(1,2,2);
camColors = lines(ncams);
hold on;
for ii = 1:ncams
    plot(1:nframes, val_reprojRMSE_cam(:,ii), 'Color', camColors(ii,:), 'LineWidth', 1.2);
end
grid on; xlabel('Frame'); ylabel('RMSE (px)');
title('Per-camera reprojection RMSE');
legend(arrayfun(@(c) sprintf('cam%d',c), camIDs, 'UniformOutput',false), 'Location','best');
xlim([1 nframes]);
sgtitle('Reprojection RMSE - CRBR MV per-frame seeding, no temporal loss');
if saveFigures
    print(gcf, fullfile(parent_folder, resultsSubdir, 'Validation_reprojRMSE_perframe.png'), '-dpng', '-r300');
end

figure('Color','w','Name','Validation: Per-body-part RMSE over time', 'Visible', onoff(showFigures), 'Position',[100 100 1200 500]);
hold on;
for j = 1:nparts
    plot(1:nframes, val_reprojRMSE_part(:,j), 'Color', partColors(j,:), 'LineWidth', 1.2);
end
grid on; xlabel('Frame'); ylabel('RMSE (px)');
title('Per-body-part reprojection RMSE over time');
legend(varnames_val, 'Location','best');
xlim([1 nframes]);
if saveFigures
    print(gcf, fullfile(parent_folder, resultsSubdir, 'Validation_reprojRMSE_perpart_overtime.png'), '-dpng', '-r300');
end
end

% Head vs body RMSE (ears count as head landmarks)
headParts = {'l_eye','r_eye','nose','l_ear','r_ear'};
bodyParts = {'spine','neck_base'};
headIdx   = find(ismember(varnames_val, headParts));
bodyIdx   = find(ismember(varnames_val, bodyParts));
rmse_head_frame = nan(nframes,1);
rmse_body_frame = nan(nframes,1);
for fidx = 1:nframes
    hVals = squeeze(val_errPix(fidx, :, headIdx)); hVals = hVals(~isnan(hVals));
    if ~isempty(hVals); rmse_head_frame(fidx) = sqrt(mean(hVals.^2)); end
    bVals = squeeze(val_errPix(fidx, :, bodyIdx)); bVals = bVals(~isnan(bVals));
    if ~isempty(bVals); rmse_body_frame(fidx) = sqrt(mean(bVals.^2)); end
end

if doFigs
figure('Color','w','Name','Validation: Head vs body RMSE', 'Visible', onoff(showFigures), 'Position',[100 100 900 420]);
hold on;
plot(1:nframes, rmse_head_frame, 'Color',[1 0 1],   'LineWidth', 1.5);
plot(1:nframes, rmse_body_frame, 'Color',[0 0.2 1], 'LineWidth', 1.5);
yline(mean(rmse_head_frame,'omitnan'), '--', 'Color',[1 0 1],   'LineWidth', 1);
yline(mean(rmse_body_frame,'omitnan'), '--', 'Color',[0 0.2 1], 'LineWidth', 1);
grid on; xlabel('Frame'); ylabel('RMSE (px)');
title('Head vs body keypoints reprojection RMSE');
legend({'Head (l\_eye, r\_eye, nose, l\_ear, r\_ear)','Body (spine, neck\_base)','Head mean','Body mean'}, 'Location','best');
xlim([1 nframes]);
if saveFigures
    print(gcf, fullfile(parent_folder, resultsSubdir, 'Validation_head_vs_body_RMSE.png'), '-dpng', '-r300');
end
end

fprintf('\nHead RMSE: mean=%.3f median=%.3f\n', mean(rmse_head_frame,'omitnan'), median(rmse_head_frame,'omitnan'));
fprintf('Body RMSE: mean=%.3f median=%.3f\n', mean(rmse_body_frame,'omitnan'), median(rmse_body_frame,'omitnan'));

% -------------------- Visual reprojection overlay (median + worst frame) --------------------
% Overlay model projections (distorted) and raw manual labels on the
% actual camera video frames for: the frame closest to median RMSE and
% the frame at maximum RMSE.

[~, repFrame]   = min(abs(val_reprojRMSE_frame - median(val_reprojRMSE_frame,'omitnan')));
repFrame        = repFrame(1);
[~, worstFrame] = max(val_reprojRMSE_frame, [], 'omitnan');
worstFrame      = worstFrame(1);

fprintf('\nOverlay figures: median LSQ %d (abs %d, RMSE %.2f px) | worst LSQ %d (abs %d, RMSE %.2f px)\n', ...
    repFrame,   startI + repFrame   - 1, val_reprojRMSE_frame(repFrame), ...
    worstFrame, startI + worstFrame - 1, val_reprojRMSE_frame(worstFrame));

ncols_fig = ceil(ncams / 2);
nrows_fig = ceil(ncams / ncols_fig);

frameSet    = [repFrame, worstFrame];
frameLabels = {'median', 'worst'};

haveOverlayVideos = doFigs && all(arrayfun(@(c) ...
    exist(fullfile(parent_folder, sprintf('cam%d_mjpeg.avi', c)),'file')==2 || ...
    exist(fullfile(parent_folder, sprintf('cam%d.mp4',        c)),'file')==2, camIDs));
if doFigs && ~haveOverlayVideos
    warning('reconPerFrame:noVideoForOverlay', ...
        'Skipping reprojection overlays: cam<N> video missing in %s', char(parent_folder));
end
if haveOverlayVideos
for fi = 1:numel(frameSet)
    fIdx   = frameSet(fi);
    fLabel = frameLabels{fi};

    [~, ~, cheetahpts_rep] = cheetah3D( ...
        xBarStore(fIdx,:), Rz, Rx, Ry, ...
        scaleHead, scaleTorsoLength);

    figure('Color','w', ...
        'Name', sprintf('Validation overlay (%s), LSQ frame %d', fLabel, fIdx), ...
        'Visible', onoff(showFigures), ...
        'Position', [60 60 420*ncols_fig 400*nrows_fig]);

    sgtitle(sprintf('Reprojection overlay - %s frame (LSQ idx %d, abs %d, RMSE %.2f px)', ...
        fLabel, fIdx, startI + fIdx - 1, val_reprojRMSE_frame(fIdx)), 'FontWeight','bold');

    for ii = 1:ncams
        cam     = camIDs(ii);
        camname = sprintf('cam%d', cam);

        videoPath = fullfile(parent_folder, sprintf('cam%d_mjpeg.avi', cam));
        if exist(videoPath, 'file') ~= 2
            videoPath = fullfile(parent_folder, sprintf('cam%d.mp4', cam));
        end
        assert(exist(videoPath, 'file') == 2, 'Video not found: %s', videoPath);

        row_offset  = find(offsets.cam == cam, 1);
        offset      = offsets.offset(row_offset);
        rawFrameIdx = startI + offset + (fIdx - 1);

        vr = VideoReader(videoPath);
        t  = max(0, min((rawFrameIdx - 1) / vr.FrameRate, vr.Duration - 1/vr.FrameRate));
        vr.CurrentTime = t;
        img = readFrame(vr);

        % Raw label table (distorted pixel coords, indexed by raw frame)
        camLabelsRaw = eval(sprintf('cam%dlabels', cam));

        subplot(nrows_fig, ncols_fig, ii);
        imshow(img); hold on;

        for j = 1:nparts
            bodypart = varnames_val{j};
            col      = partColors(j, :);

            % Undistorted coords - only for RMSE validity gating
            coord_3d    = cheetahpts_rep.(bodypart)';
            pixels_proj = world2img(coord_3d, invert(poses{ii}), intrinsicParams{ii});
            pixels_obs  = undistortedLabels.(camname).(bodypart){fIdx};

            llhIdx = find(strcmp(varsUsed, bodypart), 1);
            valid  = ~isempty(llhIdx) && ...
                     LLH(fIdx, llhIdx, ii) > 0.5 && ...
                     all(isfinite(pixels_proj)) && all(isfinite(pixels_obs)) && ...
                     all(pixels_proj > 1) && all(pixels_obs > 1);

            % Distorted coords - for plotting on the raw video frame
            pixels_proj_dist = world2img(coord_3d, invert(poses{ii}), intrinsicParams{ii}, ...
                'ApplyDistortion', true);

            xName = sprintf('%s_x', bodypart);
            yName = sprintf('%s_y', bodypart);
            if rawFrameIdx >= 1 && rawFrameIdx <= height(camLabelsRaw) && ...
                    ismember(xName, camLabelsRaw.Properties.VariableNames) && ...
                    ismember(yName, camLabelsRaw.Properties.VariableNames)
                rx = camLabelsRaw.(xName)(rawFrameIdx);
                ry = camLabelsRaw.(yName)(rawFrameIdx);
                if isfinite(rx) && isfinite(ry)
                    pixels_obs_raw = [rx ry];
                else
                    pixels_obs_raw = [NaN NaN];
                end
            else
                pixels_obs_raw = [NaN NaN];
            end

            % Always draw projection (distorted, matches raw frame)
            if all(isfinite(pixels_proj_dist)) && all(pixels_proj_dist > 1)
                plot(pixels_proj_dist(1), pixels_proj_dist(2), 'o', ...
                    'Color', col, 'MarkerSize', 8, 'LineWidth', 2);
                text(pixels_proj_dist(1) + 5, pixels_proj_dist(2), bodypart, ...
                    'Color', col, 'FontSize', 7, 'FontWeight','bold');
            end

            % Draw observation and error line only if valid
            if valid && all(isfinite(pixels_obs_raw))
                plot(pixels_obs_raw(1), pixels_obs_raw(2), '+', ...
                    'Color', col, 'MarkerSize', 10, 'LineWidth', 2);
                plot([pixels_proj_dist(1) pixels_obs_raw(1)], ...
                     [pixels_proj_dist(2) pixels_obs_raw(2)], ...
                     '-', 'Color', col, 'LineWidth', 1.5);
            end
        end

        rmse_rep = val_reprojRMSE_cam(fIdx, ii);
        title(sprintf('cam%d | raw frame %d | RMSE %.2f px', cam, rawFrameIdx, rmse_rep), ...
            'Interpreter','none');
        hold off;
    end

    if saveFigures
        savepath = fullfile(parent_folder, resultsSubdir, ...
            sprintf('Validation_overlay_%s_frame%d.png', fLabel, fIdx));
        exportgraphics(gcf, savepath, 'Resolution', 300);
    end
end
end

fprintf('\nAll outputs in: %s\n', fullfile(parent_folder, resultsSubdir));
end

% ==================== Local functions ====================
function v = getOpt(opts, name, default)
%GETOPT  Return opts.(name) if the field exists, else default.
    if isstruct(opts) && isfield(opts, name)
        v = opts.(name);
    else
        v = default;
    end
end

function s = onoff(tf)
%ONOFF  Map a logical to the 'on'/'off' char a figure 'Visible' expects.
    if tf
        s = 'on';
    else
        s = 'off';
    end
end

% Pull multiview landmarks at an absolute frame.
% seedMV    : struct with fields l_eye, r_eye, nose, l_ear, r_ear, spine, neck_base
% validMask : logical [1 x numel(landmarkNames)] in landmarkNames order
function [seedMV, validMask] = getMVseed(abs_frame, triRowOf, TRI, landmarkNames)
    seedMV    = struct();
    validMask = false(1, numel(landmarkNames));
    if ~isKey(triRowOf, double(abs_frame))
        for kk = 1:numel(landmarkNames)
            seedMV.(landmarkNames{kk}) = [NaN NaN NaN];
        end
        return
    end
    r = triRowOf(double(abs_frame));
    for kk = 1:numel(landmarkNames)
        nm = landmarkNames{kk};
        v  = TRI.pointsRaw.(nm)(r, :);
        v  = v(:).';
        seedMV.(nm) = v;
        if all(isfinite(v))
            validMask(kk) = true;
        end
    end
end

% Per-frame heading derived from multiview geometry. Returns deg2rad(180)
% if body axis (or fallback heuristics) point in -x, else 0. Validity is
% read from the seed fields directly so the landmark ordering is irrelevant.
function add = headingFromMVorPrev(seedMV, validMask, prev_yaw, hasAS, states, ~, ~) %#ok<INUSL>
    % 1) Body axis from multiview (spine -> neck_base)
    if isfield(seedMV,'spine') && isfield(seedMV,'neck_base') && ...
            all(isfinite(seedMV.spine)) && all(isfinite(seedMV.neck_base))
        bodyVec = seedMV.neck_base - seedMV.spine;
        add = double(bodyVec(1) < 0) * deg2rad(180);
        return
    end

    % 2) Head axis from multiview (centre-of-eyes -> nose)
    if isfield(seedMV,'l_eye') && isfield(seedMV,'r_eye') && isfield(seedMV,'nose') && ...
            all(isfinite(seedMV.l_eye)) && all(isfinite(seedMV.r_eye)) && all(isfinite(seedMV.nose))
        coe = (seedMV.l_eye + seedMV.r_eye) / 2;
        headVec = seedMV.nose - coe;
        add = double(headVec(1) < 0) * deg2rad(180);
        return
    end

    % 3) Previous solution yaw
    if isfinite(prev_yaw)
        add = double(abs(wrapToPi(prev_yaw)) > pi/2) * deg2rad(180);
        return
    end

    % 4) AS x_0 velocity sign as a last resort (frame 1, no prev pose)
    if hasAS && ~isempty(states)
        nS  = height(states);
        ix0 = max(1, 1 - 2);
        ix1 = min(nS, 1 + 2);
        if mean(diff(states.x_0(ix0:ix1))) < 0
            add = deg2rad(180);
        else
            add = 0;
        end
        return
    end

    add = 0;
end

% Build a 1x12 initial guess and a "filled" seed struct.
% Missing landmarks (validMask == false) are replaced with the previous
% frame's kinematic-model points so ASseedResidualFromModel works as-is.
%
% HEAD orientation is seeded by a rotation-only Kabsch/Umeyama fit of the
% template head landmarks (eyes, nose, ears) to the triangulated head points,
% giving a CONTINUOUS full-360 deg head-yaw seed that is independent of the
% body axis. The old binary 0/180 deg `addition` placed the head seed in the
% wrong yaw half-circle during the turn, where the head and body point in very
% different directions; the eyes+nose+ear residuals then trapped the per-frame
% solver in the flipped basin. `addition_i` is now used for the head only as a
% fallback when fewer than 3 head landmarks are valid.
%
% TORSO orientation is seeded continuously from the triangulated spine->neck
% axis (a directed vector, so full-360 deg with no +/-180 ambiguity), with roll
% held nominal; `addition_i` is used for the torso only as a fallback when spine
% or neck_base is not validly triangulated this frame. The old binary 0/180 deg
% `addition` could land the torso seed in the wrong yaw half-circle during a
% turn (the per-frame solve then settled ~180 deg off).
function [params0_i, pHead0_i, torso0_i, sf] = buildParams0( ...
        seedMV, validMask, prev_pts, addition_i, scaleHead, scaleTorsoLength, Rx, Ry, Rz)

    names = {'l_eye','r_eye','nose','l_ear','r_ear','spine','neck_base'};
    sf = struct();
    for k = 1:numel(names)
        nm = names{k};
        if validMask(k)
            sf.(nm) = seedMV.(nm);
        else
            sf.(nm) = prev_pts.(nm)(:).';
        end
    end

    % ---- HEAD: continuous orientation via rotation-only Kabsch fit ----
    % Template head landmarks in head-local coords (must match cheetah3D).
    headNames = {'l_eye','r_eye','nose','l_ear','r_ear'};
    tmplHead  = scaleHead * [ ...
         0       0.03    0      ; ...   % l_eye
         0      -0.03    0      ; ...   % r_eye
         0.055   0      -0.055  ; ...   % nose
        -0.0725  0.0608  0.0441 ; ...   % l_ear
        -0.0725 -0.0608  0.0441 ];      % r_ear

    headValid = validMask(1:5);          % validMask is in `names` order
    P  = tmplHead(headValid, :);
    vi = find(headValid);
    Q  = zeros(numel(vi), 3);
    for j = 1:numel(vi)
        Q(j,:) = sf.(headNames{vi(j)})(:).';
    end

    if size(P,1) >= 3 && rank(P - mean(P,1)) >= 2
        [R_head, t_head] = kabschR(P, Q);
        [hx, hy, hz]     = eulZXY(R_head);
        pHead0_i         = t_head(:);
    else
        % Sparse head coverage: fall back to the binary-heading seed.
        R_init   = Rz(addition_i) * Rx(0) * Ry(0);
        or2nose  = R_init * [scaleHead*0.055; 0; -scaleHead*0.055];
        pHead0_i = sf.nose(:) - or2nose;
        hx = 0; hy = 0; hz = addition_i;
    end

    % ---- TORSO: continuous orientation from the triangulated body axis ----
    % The torso has only two landmarks (spine, neck_base), so a full Kabsch fit
    % cannot determine a rotation (and must not: roll about the spine->neck axis
    % is a degenerate DOF that does not affect the body forward axis / pitch /
    % yaw the analysis uses). Instead build a zero-roll "level" rotation that
    % maps local +x onto the triangulated spine->neck vector, then extract Euler
    % angles with the SAME eulZXY helper the head uses (convention-consistent
    % with cheetah3D's R1_3 = Rz*Rx*Ry). Because the body vector is a
    % DIRECTED 3-D vector (neck minus spine), the heading is continuous over the
    % full 360 deg with NO +/-180 ambiguity -- this resolves the branch that the
    % old binary sign(neck_base.x - spine.x) seed flipped during turns.
    %
    % Requires both spine and neck_base to be validly triangulated this frame
    % (validMask is in `names` order: spine = 6, neck_base = 7). When they are
    % not, fall back to the original binary-heading (`addition_i`) seed.
    bodyVec   = sf.neck_base(:) - sf.spine(:);
    side      = cross([0;0;1], bodyVec);   % horizontal local-y => roll ~ 0 (level)
    torsoCont = validMask(6) && validMask(7) && ...
                norm(bodyVec) > 1e-6 && norm(side) > 1e-6;

    if torsoCont
        u   = bodyVec / norm(bodyVec);     % forward (spine -> neck_base)
        sde = side / norm(side);           % local-y image (horizontal)
        upl = cross(u, sde);               % local-z image
        R_torso = [u, sde, upl];           % columns = images of local x,y,z
        [tx2, ty2, tz2] = eulZXY(R_torso);
    else
        % Sparse / degenerate body axis: fall back to the binary-heading seed.
        R_torso = Rz(addition_i) * Rx(0) * Ry(0);
        tx2 = 0; ty2 = 0; tz2 = addition_i;
    end

    or2midspine_world = R_torso * [-scaleTorsoLength*0.33; 0; scaleTorsoLength*0.10];
    torso0_i          = sf.spine(:) - or2midspine_world;

    params0_i        = zeros(1,12);
    params0_i(1:3)   = [hx hy hz];
    params0_i(4:6)   = pHead0_i.';
    params0_i(7:9)   = [tx2 ty2 tz2];
    params0_i(10:12) = torso0_i.';
end

% Rotation-only Kabsch/Umeyama fit. Returns R (local->world) and t such that
% Q(i,:).' ~= R*P(i,:).' + t for matched rows of P (head-local template) and
% Q (triangulated world points). Mirrors the SVD alignment used in
% templateGeometryFromTriangulation.m.
function [R, t] = kabschR(P, Q)
    mP = mean(P, 1);
    mQ = mean(Q, 1);
    Pc = P - mP;
    Qc = Q - mQ;
    H  = Pc.' * Qc;                       % 3x3 cross-covariance
    [U, ~, V] = svd(H);
    d = sign(det(V * U.'));
    if d == 0, d = 1; end
    R = V * diag([1 1 d]) * U.';          % proper rotation (det = +1)
    t = mQ.' - R * mP.';
end

% Euler-angle extraction for the script's convention R = Rz(tz)*Rx(tx)*Ry(ty),
% where Rx, Rz are the transposed (negated-angle) forms and Ry is standard (see
% Rx/Ry/Rz in the main script). This is the inverse of cheetah3D's R1_2.
function [tx, ty, tz] = eulZXY(R)
    tx = atan2(-R(3,2), sqrt(R(1,2)^2 + R(2,2)^2));
    tz = atan2( R(1,2), R(2,2));
    ty = atan2(-R(3,1), R(3,3));
end

% Raw (non-robust) sigma-normalised reprojection RMSE for one solution, used to
% pick the better multi-start result. robustOpt is omitted so cheetahError
% uses plain residuals; useTemporal=false and params_before is unused.
function rmse = reprojRMSEraw(p, Rz, Rx, Ry, scaleHead, scaleTorsoLength, ...
        undistLabels, poses, intrinsicParams, LLH, frameNum, variances, camIDs)
    E = cheetahError(p, Rz, Rx, Ry, scaleHead, scaleTorsoLength, ...
        undistLabels, poses, intrinsicParams, LLH, frameNum, variances, p, camIDs, false);
    ncams = numel(camIDs);
    R_all = [];
    for c = 1:ncams
        R_all = [R_all; E(:,:,c)];
    end
    rmse = sqrt(mean(R_all(:).^2));
end

% Draw the cheetah head/body skeleton (eye line, centre-of-eyes -> nose, and
% spine -> neck_base) into the current 3D axes. Segments whose endpoints are
% missing or non-finite are skipped, so a partial pose such as AS (no ears) and
% any landmark gaps are handled gracefully.
function drawSkeleton3D(pts, color, lineStyle, lw)
    seg = @(p, q) plot3([p(1) q(1)], [p(2) q(2)], [p(3) q(3)], ...
        'LineStyle', lineStyle, 'Color', color, 'LineWidth', lw);

    okEyes = isfield(pts,'l_eye') && isfield(pts,'r_eye') && ...
             all(isfinite(pts.l_eye(:))) && all(isfinite(pts.r_eye(:)));
    if okEyes
        le = pts.l_eye(:).'; re = pts.r_eye(:).';
        seg(le, re);
        if isfield(pts,'nose') && all(isfinite(pts.nose(:)))
            coe = (le + re) / 2;
            seg(coe, pts.nose(:).');
        end
    end
    if isfield(pts,'spine') && isfield(pts,'neck_base') && ...
            all(isfinite(pts.spine(:))) && all(isfinite(pts.neck_base(:)))
        seg(pts.spine(:).', pts.neck_base(:).');
    end
end
