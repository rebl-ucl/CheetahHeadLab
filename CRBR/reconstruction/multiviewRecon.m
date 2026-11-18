function multiviewRecon(parent_folder, camIDs, opts)
% multiviewRecon  CRBR multiview triangulation from manual labels.
%
%   multiviewRecon(parent_folder, camIDs, opts)
%
% Uses all visible cameras for each body part and frame, undistorted points, a
% camera-pose table + pointTrack + triangulateMultiview. No kinematic model /
% lsqnonlin / temporal regulariser. A leave-one-out per-camera trust pass
% (cameraTrustLOO) rejects consensus-outlier cameras per frame.
%
% Writes <parent_folder>/triangulatedPointsRaw_multiview.mat (the seed for
% reconPerFrame) plus optional med/lpf variants and a camera-rejection log.
% triangulateMultiview ignores lens distortion, so points are undistorted first;
% the cameraPoses ViewId must match the pointTrack viewIDs exactly.
%
% opts fields (all optional; defaults are headless + fast for batch):
%   showFigures(false) createVideo(false) createOrthoVideo(false)
%   plotRawResults/plotMedResults/plotLPFResults(false) pointSource("raw")
%   applyMedianFilter/applyLPF(false) medianWindow([]) fc(5) fs(90)
%   inspectFrame(10) inspectCam(2) doInspect(false) trustOpt(struct)
%
% Called with no arguments it runs a default dataset (standalone debugging).

if nargin < 1 || isempty(parent_folder)
    parent_folder = fullfile(pwd, 'data', 'pitch_tracking', 'phantomRun_20190307');
end
if nargin < 2 || isempty(camIDs)
    camIDs = [1 2 3 4 5 6];
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

ncams = numel(camIDs);

% Reconstruct ears (not shoulders) so the result can seed the full-head-
% trajectory estimation. Order MUST match varsUsed in reconPerFrame.
varsUsed = {'l_eye', 'r_eye', 'nose', 'l_ear', 'r_ear', 'spine', 'neck_base'};
nBody = numel(varsUsed);

visibilityThreshold = 0.5;   % LLH > threshold => visible
minViewsToTriangulate = 2;

% Leave-one-out per-camera trust (cameraTrustLOO.m). A camera is REJECTED for
% a frame when its leave-one-out reprojection residual is a consensus outlier
% (> C_rel x median of peers AND > floorPx), sustained for H frames.
trustOpt = getOpt(opts, 'trustOpt', struct( ...
    'C_rel',          3, ...   % reject if LOO residual > C_rel x median(peers)
    'floorPx',        8, ...   % ... and > this many px (keeps all cams on good frames)
    'minPresentCams', 4, ...   % only reject when >= this many cams present
    'keepAtLeast',    3, ...   % never drop the kept set below this many cams
    'H',              2, ...   % hysteresis: flip kept<->rejected after H frames
    'minCamsLOO',     3));      % a landmark needs >= this many views for a LOO test

% Figure / video show/save convention (shared across the pipeline):
%   showFigures - pop up diagnostic figure windows
%   saveFigures - write the diagnostic report PNGs to reconResults_triangulation/
%   showVideos  - render the reconstruction video live while writing
%   saveVideos  - write the reconstruction video (createVideo) to .avi
showFigures      = getOpt(opts, 'showFigures',       false);
saveFigures      = getOpt(opts, 'saveFigures',       true);
showVideos       = getOpt(opts, 'showVideos',        false);
saveVideos       = getOpt(opts, 'saveVideos',        false);
createVideo      = getOpt(opts, 'createVideo',       saveVideos);   % primary 3D recon video
createOrthoVideo = getOpt(opts, 'createOrthoVideo',  false);        % secondary ortho video (opt-in)
% Build the diagnostic figures whenever they will be shown or saved; each
% figure's on-screen visibility tracks showFigures, the PNG write tracks
% saveFigures (so unattended batch runs save without popping windows).
doFigs           = showFigures || saveFigures;

plotRawResults   = getOpt(opts, 'plotRawResults',    false);
plotMedResults   = getOpt(opts, 'plotMedResults',    false);
plotLPFResults   = getOpt(opts, 'plotLPFResults',    false);

pointSource      = getOpt(opts, 'pointSource',       "raw");   % "raw", "med", or "lpf"

applyMedianFilter= getOpt(opts, 'applyMedianFilter', false);
applyLPF         = getOpt(opts, 'applyLPF',          false);

medianWindow     = getOpt(opts, 'medianWindow',      []);   % auto if empty
fc               = getOpt(opts, 'fc', 5);               % Hz
fs               = getOpt(opts, 'fs', 90);              % Hz

inspectFrame     = getOpt(opts, 'inspectFrame', 10);   % index inside reconstruction range
inspectCam       = getOpt(opts, 'inspectCam',   2);    % actual camera number in camIDs
doInspect        = getOpt(opts, 'doInspect',    false);
%%
%[text] ## -------------------- Output folder --------------------
reconFolder = fullfile(parent_folder, "reconResults_triangulation");
if ~exist(reconFolder, 'dir')
    mkdir(reconFolder);
end
%%
%[text] ## -------------------- Load labels + LLH --------------------
labelTables = cell(1, ncams);
llhTables   = cell(1, ncams);

for k = 1:ncams  
    camNum = camIDs(k);

    labelFile = fullfile(parent_folder, sprintf("cam%d_labels_wide.csv", camNum));
    llhFile   = fullfile(parent_folder, sprintf("cam%d_llh_wide.csv", camNum));

    assert(exist(labelFile,'file')==2, "Missing file: %s", labelFile);
    assert(exist(llhFile,'file')==2,   "Missing file: %s", llhFile);

    labelTables{k} = readtable(labelFile);
    llhTables{k}   = readtable(llhFile);

    fprintf("cam%d: labels size = %s, llh size = %s\n", ...  
        camNum, mat2str(size(labelTables{k})), mat2str(size(llhTables{k})));  
end  
%%
%[text] ## -------------------- Offsets --------------------
osfile = dir(fullfile(parent_folder, "AcinoSet_data", "*offset*"));
if ~isempty(osfile)
    offsetsTbl = readtable(fullfile(parent_folder, "AcinoSet_data", osfile(1).name), ...
        VariableNamingRule="preserve");
else
    offsetsTbl = table(zeros(ncams,1), 'VariableNames', {'offset'});
end

assert(ismember('offset', offsetsTbl.Properties.VariableNames), ...
    'Offsets table must contain column "offset".');

offsets = offsetsTbl.offset;

for k = 1:ncams  
    fprintf('offset(cam%d) = %g\n', camIDs(k), offsets(k));   
end  
%[text] 
%[text] ## -------------------- Load camera parameters --------------------
jsonCandidates = dir(fullfile(parent_folder, "AcinoSet_data", "*cam_scene_sba.json"));
assert(~isempty(jsonCandidates), "No *cam_scene_sba.json found.");

wantedName = sprintf("%d_cam_scene_sba.json", ncams);
match = find(strcmp({jsonCandidates.name}, wantedName), 1);

assert(~isempty(match), "Could not find %s in AcinoSet_data.", wantedName);

extrinsicsFile = fileread(fullfile(parent_folder, "AcinoSet_data", jsonCandidates(match).name));

extrinsics = jsondecode(extrinsicsFile);

imageSize = [extrinsics.camera_resolution(2), extrinsics.camera_resolution(1)];

intrinsicParams = cell(1, ncams);   % keep as returned by cameraIntrinsicsFromOpenCV
posesPlot       = cell(1, ncams);   % camera->world for plotting
%%
%[text] ## -------------------- Camera pose table for triangulateMultiview --------------------
cameraPoseTable = table( ...
    zeros(ncams,1,'uint32'), ...
    repmat(rigidtform3d, ncams, 1), ...
    'VariableNames', {'ViewId','AbsolutePose'});

for k = 1:ncams
    c = extrinsics.cameras(k);

    R_cw = c.r;       % world -> camera
    t_cw = c.t(:);    % 3x1

    intrin = cameraIntrinsicsFromOpenCV(c.k, c.d, imageSize);
    intrinsicParams{k} = intrin;

    % camera -> world pose
    R_wc = R_cw.';
    C_w  = -(R_cw.' * t_cw);
    T_wc = rigidtform3d(R_wc, C_w.');

    posesPlot{k} = T_wc;

    % IMPORTANT: exact variable name is ViewId
    cameraPoseTable.ViewId(k) = uint32(camIDs(k));
    cameraPoseTable.AbsolutePose(k) = T_wc;

    % projection matrix for reprojection diagnostics
    % K = intrin.K';
    % P_all{k} = K * [R_cw, t_cw];
end

cameraPoses = cameraPoseTable;

if showFigures
    figure('Color','w');
    plotCameraConfigLT6(posesPlot, 0.1, camIDs);
    title('Camera configuration');
    axis vis3d
    grid on
end
%%
%[text] ## -------------------- Reconstruction frame range --------------------
refLabels = labelTables{1};

limitfile = dir(fullfile(parent_folder, "*limit*"));
if ~isempty(limitfile) 
    lims = readtable(fullfile(parent_folder, limitfile(1).name), ...
        VariableNamingRule="preserve");
    startI = lims{1,1};
    endI   = lims{1,2};
    fprintf("Using limit file frame range: %d -> %d\n", startI, endI); 
else
    startI = 1;
    endI   = height(refLabels);
    fprintf("No limit file found. Using all available frames: %d -> %d\n", startI, endI);
end 

startI = max(startI, 1);
endI   = min(endI, height(refLabels));
assert(startI <= endI, "Invalid reconstruction frame range.");

nframes   = endI - startI + 1;
absFrames = (startI:endI).';

fprintf("Final frame range: %d -> %d (nframes=%d)\n", startI, endI, nframes); 
%%
%[text] ## -------------------- Crop + undistort labels --------------------
undistortedLabels = struct();
rawLabelsCropped  = struct();
likelihoods       = struct();

for k = 1:ncams
    camNum  = camIDs(k);
    camname = sprintf('cam%d', camNum);
    intrin  = intrinsicParams{k};
    offset  = offsets(k);

    labelsTbl = labelTables{k};
    llhTbl    = llhTables{k};

    if ismember("frame", labelsTbl.Properties.VariableNames)
        labelsTbl = removevars(labelsTbl, "frame");
    end
    if ismember("frame", llhTbl.Properties.VariableNames)
        llhTbl = removevars(llhTbl, "frame");
    end

    rowStart = startI + offset;
    rowEnd   = endI   + offset;

    assert(rowStart >= 1, "cam%d: rowStart=%d < 1", camNum, rowStart);
    assert(rowEnd <= height(labelsTbl), ...
        "cam%d: rowEnd=%d exceeds labels rows=%d", camNum, rowEnd, height(labelsTbl));
    assert(rowEnd <= height(llhTbl), ...
        "cam%d: rowEnd=%d exceeds llh rows=%d", camNum, rowEnd, height(llhTbl));

    labelsCrop = labelsTbl(rowStart:rowEnd, :);
    llhCrop    = llhTbl(rowStart:rowEnd, :);

    labelsKeep = table();
    llhKeep    = table();

    for b = 1:nBody
        bp = varsUsed{b};

        xName = sprintf('%s_x', bp);
        yName = sprintf('%s_y', bp);
        lName = sprintf('%s_LLH', bp);

        assert(ismember(xName, labelsCrop.Properties.VariableNames), ...
            "cam%d missing %s", camNum, xName);
        assert(ismember(yName, labelsCrop.Properties.VariableNames), ...
            "cam%d missing %s", camNum, yName);
        assert(ismember(lName, llhCrop.Properties.VariableNames), ...
            "cam%d missing %s", camNum, lName);

        labelsKeep.(xName) = labelsCrop.(xName);
        labelsKeep.(yName) = labelsCrop.(yName);
        llhKeep.(bp)       = llhCrop.(lName);
    end

    labelsUD = table();
    for b = 1:nBody
        bp = varsUsed{b};
        xName = sprintf('%s_x', bp);
        yName = sprintf('%s_y', bp);

        pts = [labelsKeep.(xName), labelsKeep.(yName)];
        valid = all(isfinite(pts),2);

        ptsUD = nan(size(pts));
        if any(valid)
            ptsUD(valid,:) = undistortPoints(pts(valid,:), intrin);
        end

        labelsUD.(xName) = ptsUD(:,1);
        labelsUD.(yName) = ptsUD(:,2);
    end

    rawLabelsCropped.(camname)  = labelsKeep;
    undistortedLabels.(camname) = labelsUD;
    likelihoods.(camname)       = llhKeep;
end
%%
%[text] ## -------------------- Intrinsics array for triangulateMultiview --------------------
%[text] Many MATLAB examples use one shared intrinsics object when all views come from the same physical camera model. Here each camera can have its own. We pass a vector / cell-expanded array in the same order as cameraPoses.
%[text] If your MATLAB release rejects the object type here, this is where the compatibility issue will surface.
intrinsicsForTriangulation = [intrinsicParams{:}];
%%
%[text] ## -------------------- Multiview triangulation --------------------
points3D_raw      = nan(nframes, nBody, 3);   % FINAL: triangulated from TRUSTED cams only
points3D_allcam   = nan(nframes, nBody, 3);   % diagnostic: all visible cams (pre-rejection)
triangOK          = false(nframes, nBody);
numViewsUsed      = zeros(nframes, nBody);     % # visible views (before trust rejection)
numViewsKept      = zeros(nframes, nBody);     % # trusted views used for points3D_raw
reprojRMSE_body   = nan(nframes, nBody);       % from triangulateMultiview output (trusted cams)
validIndex_body   = false(nframes, nBody);

% custom reprojection diagnostics across ALL visible cams (raw, pre-rejection)
reprojErrPerObs         = nan(nframes, nBody, ncams, 2);  % raw dx, dy in undistorted px

reprojRMSE_frame        = nan(nframes,1);   % raw px RMSE (all visible cams)
coverage_pairs          = zeros(nframes,1);

% per-(frame,bodypart) observation cache for the leave-one-out camera-trust pass
obsCell = cell(nframes, nBody);

% ---- Pass 1: all-camera triangulation (diagnostics) + capture observations ----
for f = 1:nframes
    for b = 1:nBody
        bp = varsUsed{b};
        xName = sprintf('%s_x', bp);
        yName = sprintf('%s_y', bp);

        viewIDs = zeros(0,1,'uint32');
        pts2D   = zeros(0,2);

        for k = 1:ncams
            camNum  = camIDs(k);
            camname = sprintf('cam%d', camNum);

            uv = [undistortedLabels.(camname).(xName)(f), ...
                  undistortedLabels.(camname).(yName)(f)];
            llh = likelihoods.(camname).(bp)(f);

            if all(isfinite(uv)) && isfinite(llh) && llh > visibilityThreshold
                viewIDs(end+1,1) = uint32(camNum);
                pts2D(end+1,:)   = uv;
            end
        end

        numViewsUsed(f,b) = numel(viewIDs);

        % cache observations for the camera-trust (LOO) pass below
        obsCell{f,b} = struct('viewIDs', viewIDs, 'pts2D', pts2D);

        if numel(viewIDs) < minViewsToTriangulate
            continue
        end

        % one point observed in multiple views -> one pointTrack object
        track = pointTrack(viewIDs.', pts2D);

        try
            [X, ~, validIdx] = triangulateMultiview(track, cameraPoses, intrinsicsForTriangulation);
        catch ME
            warning('triangulateMultiview failed at frame %d bodypart %s: %s', ...
                f, bp, ME.message);
            continue
        end

        if isempty(X) || ~all(isfinite(X)) || isempty(validIdx) || ~validIdx(1)
            continue
        end

        points3D_allcam(f,b,:) = X(1,:);

        % ---- custom reprojection diagnostics in undistorted coordinates (all cams) ----
        Xw = X(1,:);
        for k = 1:ncams
            camNum = camIDs(k);
            if ismember(camNum, viewIDs)
                % Reproject to undistorted image coords (matches undistorted pts2D)
                uv_proj = world2img(Xw, invert(posesPlot{k}), intrinsicParams{k});
                uv_obs  = pts2D(viewIDs == camNum, :);
                reprojErrPerObs(f,b,k,:) = uv_proj - uv_obs;
            end
        end
    end

    E = squeeze(reprojErrPerObs(f,:,:,:));        % raw px residuals (all cams)
    E = reshape(E, [], 2);
    validRows = all(isfinite(E),2);

    if any(validRows)
        coverage_pairs(f)   = sum(validRows);
        reprojRMSE_frame(f) = sqrt(mean(E(validRows,:).^2, 'all'));
    end
end

% ---- Pass 2: leave-one-out camera-trust decision (per frame, per camera) ----
% A close camera carrying a coherent calibration / sub-frame-sync bias keeps a
% LOW in-sample residual (the joint fit moves toward it), so it must be judged by
% leave-one-out: triangulate each landmark from the OTHER cams, reproject into the
% held-out cam, and flag cameras whose LOO residual is a consensus outlier.
[camKeptMask, camLOOResid] = cameraTrustLOO( ...
    obsCell, cameraPoses, intrinsicsForTriangulation, posesPlot, intrinsicParams, ...
    camIDs, trustOpt);

% ---- Pass 3: re-triangulate each landmark from TRUSTED cameras only ----
for f = 1:nframes
    keptCamNums = camIDs(camKeptMask(f,:));
    for b = 1:nBody
        o = obsCell{f,b};
        if isempty(o), continue; end

        keepRows = ismember(double(o.viewIDs(:)), keptCamNums);
        vID = o.viewIDs(keepRows);
        p2D = o.pts2D(keepRows,:);
        numViewsKept(f,b) = numel(vID);

        if numel(vID) < minViewsToTriangulate
            continue
        end

        track = pointTrack(vID.', p2D);
        try
            [X, reprojErr, validIdx] = triangulateMultiview(track, cameraPoses, intrinsicsForTriangulation);
        catch ME
            warning('trusted-cam triangulateMultiview failed at frame %d bodypart %s: %s', ...
                f, varsUsed{b}, ME.message);
            continue
        end
        if isempty(X) || ~all(isfinite(X)) || isempty(validIdx) || ~validIdx(1)
            continue
        end

        points3D_raw(f,b,:) = X(1,:);
        triangOK(f,b)       = true;
        validIndex_body(f,b)= validIdx(1);
        reprojRMSE_body(f,b)= reprojErr(1);
    end
end

% Per-camera, per-frame in-sample reprojection RMSE (from all-cam reprojErrPerObs),
% used for the transparency log/plot alongside the LOO score and the kept mask.
reprojRMSE_cam_perframe = nan(nframes, ncams);
for f = 1:nframes
    for k = 1:ncams
        Ec    = squeeze(reprojErrPerObs(f,:,k,:));  % [nBody x 2]
        valid = all(isfinite(Ec), 2);
        if any(valid)
            reprojRMSE_cam_perframe(f,k) = sqrt(mean(Ec(valid,:).^2, 'all'));
        end
    end
end
%%
%[text] ## -------------------- Summary --------------------
fprintf('\n=== Multiview triangulation summary ===\n'); 
fprintf('Frames: %d\n', nframes); 
fprintf('Body parts: %d\n', nBody); 
fprintf('Mean triangulated body parts / frame: %.2f / %d\n', ...  
    mean(sum(triangOK,2)), nBody);  
fprintf('Mean visible observations / frame: %.2f\n', mean(coverage_pairs)); 

validRMSE = reprojRMSE_frame(isfinite(reprojRMSE_frame));
if ~isempty(validRMSE) 
    fprintf('Custom reprojection RMSE (undistorted px): mean=%.4f, median=%.4f, p95=%.4f, max=%.4f\n', ... 
        mean(validRMSE), median(validRMSE), prctile(validRMSE,95), max(validRMSE)); 
end 

validBodyErr = reprojRMSE_body(isfinite(reprojRMSE_body));
if ~isempty(validBodyErr) 
    fprintf('triangulateMultiview returned reprojection error: mean=%.4f, median=%.4f, p95=%.4f, max=%.4f\n', ... 
        mean(validBodyErr), median(validBodyErr), prctile(validBodyErr,95), max(validBodyErr)); 
end 
%%
%[text] ## -------------------- Save raw result --------------------
pointsRaw = packPointsStruct(points3D_raw, varsUsed);

save(fullfile(parent_folder, "triangulatedPointsRaw_multiview.mat"), ...
    "points3D_raw", "pointsRaw", "varsUsed", "camIDs", "absFrames", ...
    "triangOK", "numViewsUsed", "numViewsKept", "validIndex_body", ...
    "reprojErrPerObs", ...
    "reprojRMSE_body", ...
    "reprojRMSE_frame", ...
    "reprojRMSE_cam_perframe", ...
    "coverage_pairs", ...
    "points3D_allcam", ...
    "camKeptMask", "camLOOResid", "trustOpt");
%%
%[text] ## -------------------- Camera-trust transparency log --------------------
%[text] Per-frame record of which cameras were kept vs rejected by the leave-one-out
%[text] trust pass, with each camera's LOO residual and in-sample reprojection RMSE,
%[text] so the rejection process is fully auditable.
nKeptPerFrame = sum(camKeptMask, 2);

trustTbl = table((1:nframes).', absFrames(:), nKeptPerFrame, ...
    'VariableNames', {'frame','absFrame','nKept'});
for k = 1:ncams
    cid = camIDs(k);
    trustTbl.(sprintf('cam%d_LOOresid',   cid)) = camLOOResid(:,k);
    trustTbl.(sprintf('cam%d_reprojRMSE', cid)) = reprojRMSE_cam_perframe(:,k);
    trustTbl.(sprintf('cam%d_kept',       cid)) = double(camKeptMask(:,k));
end
trustCsv = fullfile(reconFolder, "camRejectionLog_multiview.csv");
writetable(trustTbl, trustCsv);

fprintf('\n=== Camera-trust (leave-one-out) summary ===\n');
fprintf('Per-frame camera rejection log: %s\n', trustCsv);
for k = 1:ncams
    cid        = camIDs(k);
    rejected_k = sum(~camKeptMask(:,k));
    medLOO     = median(camLOOResid(isfinite(camLOOResid(:,k)),k));
    fprintf('  cam%d: rejected in %d / %d frames (median LOO residual %.2f px)\n', ...
        cid, rejected_k, nframes, medLOO);
end
rejFrames = find(nKeptPerFrame < ncams);
if isempty(rejFrames)
    fprintf('  No cameras rejected on any frame (all consensus-consistent).\n');
else
    fprintf('  Frames with >=1 rejected camera: %d (abs %s)\n', ...
        numel(rejFrames), mat2str(absFrames(rejFrames).'));
end

if doFigs
% Per-camera reprojection RMSE over time, with LOO rejection markers
figure('Color','w','Name','Per-camera reproj RMSE + LOO rejections', 'Visible', onoff(showFigures));
hold on;
trustCamColors = lines(ncams);
for k = 1:ncams
    plot(1:nframes, reprojRMSE_cam_perframe(:,k), '-', ...
        'Color', trustCamColors(k,:), 'LineWidth', 1.2);
end
for k = 1:ncams
    rej = find(~camKeptMask(:,k));
    if ~isempty(rej)
        plot(rej, reprojRMSE_cam_perframe(rej,k), 'x', ...
            'Color', trustCamColors(k,:), 'MarkerSize', 9, 'LineWidth', 1.6);
    end
end
grid on; xlabel('Frame index'); ylabel('In-sample reproj RMSE (px)');
title('Per-camera reprojection RMSE (x = rejected by leave-one-out trust)');
legend(arrayfun(@(c) sprintf('cam%d', c), camIDs, 'UniformOutput', false), 'Location','best');
if saveFigures
    exportgraphics(gcf, fullfile(reconFolder, "Per-camera_reproj_RMSE_LOO_rejections.png"), 'Resolution', 300);
end
%%
%[text] ## -------------------- Diagnostics --------------------
figure('Color','w', 'Visible', onoff(showFigures));
plot(reprojRMSE_frame, 'LineWidth', 1.2);
grid on
xlabel('Frame index');
ylabel('RMSE (undistorted px)');
title('Per-frame reprojection RMSE (multiview)');
if saveFigures
    exportgraphics(gcf, fullfile(reconFolder, "Per-frame_reprojection_RMSE_multiview.png"), 'Resolution', 300);
end

figure('Color','w', 'Visible', onoff(showFigures));
plot(coverage_pairs, 'LineWidth', 1.2);
grid on
xlabel('Frame index');
ylabel('# reprojection observation pairs');
title('Per-frame observation coverage');
if saveFigures
    exportgraphics(gcf, fullfile(reconFolder, "Per-frame_observation_coverage_multiview.png"), 'Resolution', 300);
end

figure('Color','w', 'Visible', onoff(showFigures));
imagesc(numViewsUsed');
colorbar
xlabel('Frame index');
ylabel('Body part');
yticks(1:nBody);
yticklabels(varsUsed);
title('Number of views used per body part');
if saveFigures
    exportgraphics(gcf, fullfile(reconFolder, "Num_views_used_per_bodypart_multiview.png"), 'Resolution', 300);
end
end
%%
%[text] ## -------------------- Filtering --------------------
points3D_med = points3D_raw;
points3D_lpf = points3D_raw;

if isempty(medianWindow)
    medianWindow = round(nframes * 0.0625);
    medianWindow = max(medianWindow, 3);
    if mod(medianWindow,2)==0
        medianWindow = medianWindow + 1;
    end
end

if applyMedianFilter
    points3D_med = filter3DTrajectories(points3D_raw, 'median', medianWindow, []);
    pointsMed = packPointsStruct(points3D_med, varsUsed);
    save(fullfile(parent_folder, "triangulatedPointsMed_multiview.mat"), ...
        "points3D_med", "pointsMed", "varsUsed", "camIDs", "absFrames", "medianWindow");
end

if applyLPF
    points3D_lpf = filter3DTrajectories(points3D_raw, 'lpf', [], struct('fc',fc,'fs',fs));
    pointsLPF = packPointsStruct(points3D_lpf, varsUsed);
    save(fullfile(parent_folder, "triangulatedPointsLPF_multiview.mat"), ...
        "points3D_lpf", "pointsLPF", "varsUsed", "camIDs", "absFrames", "fc", "fs");
end
%%
%[text] ## -------------------- Trajectory plots --------------------
if plotRawResults
    plotPointTrajectories(points3D_raw, varsUsed, absFrames, 'Raw triangulated 3D trajectories (multiview)');
    exportgraphics(gcf, fullfile(reconFolder, "Triangulated_3D_trajectories_raw_multiview.png"), 'Resolution', 300);
end

if plotMedResults
    plotPointTrajectories(points3D_med, varsUsed, absFrames, 'Median-filtered triangulated 3D trajectories (multiview)');
    exportgraphics(gcf, fullfile(reconFolder, "Triangulated_3D_trajectories_med_multiview.png"), 'Resolution', 300);
end

if plotLPFResults
    plotPointTrajectories(points3D_lpf, varsUsed, absFrames, 'Low-pass filtered triangulated 3D trajectories (multiview)');
    exportgraphics(gcf, fullfile(reconFolder, "Triangulated_3D_trajectories_lpf_multiview.png"), 'Resolution', 300);
end
%%
%[text] ## -------------------- Select source for video --------------------
switch lower(pointSource)
    case "raw"
        points3D_video = points3D_raw;
        videoTag = "raw";
    case "med"
        points3D_video = points3D_med;
        videoTag = "med";
    case "lpf"
        points3D_video = points3D_lpf;
        videoTag = "lpf";
    otherwise
        error('pointSource must be "raw", "med", or "lpf"');
end
%%
%[text] ## -------------------- 3D video --------------------
if createVideo
    outName = fullfile(reconFolder, "TriangulationRecon_multiview_" + videoTag + ".avi");
    vw = VideoWriter(outName, 'Motion JPEG AVI');
    vw.Quality = 100;
    vw.FrameRate = 20;
    open(vw);

    fg = figure('Color','w','Visible',onoff(showVideos));
    set(fg, 'Position', [80 80 1600 900]);

    for f = 1:nframes
        clf
        plotCameraConfigLT6(posesPlot, 0.1, camIDs); hold on;

        P = squeeze(points3D_video(f,:,:)); % [nBody x 3]
        valid = all(isfinite(P),2);

        if any(valid)
            plotTriangulatedCheetah(P, varsUsed);
        end

        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(sprintf('Multiview triangulated 3D reconstruction (%s), absolute frame %d', ...
            videoTag, absFrames(f)));
        axis equal
        grid on
        view(3)
        drawnow

        writeVideo(vw, getframe(fg));
    end

    close(vw);
    close(fg);

    fprintf('3D video saved: %s\n', outName);
end
%%
%[text] ## -------------------- Orthographic video --------------------
if createOrthoVideo
    outName = fullfile(reconFolder, "TriangulationRecon_ortho_multiview_" + videoTag + ".avi");
    vw = VideoWriter(outName, 'Motion JPEG AVI');
    vw.Quality = 100;
    vw.FrameRate = 20;
    open(vw);

    fg = figure('Color','w','Visible',onoff(showVideos));
    set(fg, 'Position', [80 80 1800 600]);

    for f = 1:nframes
        clf
        P = squeeze(points3D_video(f,:,:)); % [nBody x 3]
        valid = all(isfinite(P),2);

        if any(valid)
            xs = P(valid,1); ys = P(valid,2); zs = P(valid,3);
            xlims = [min(xs)-0.08, max(xs)+0.08];
            ylims = [min(ys)-0.08, max(ys)+0.08];
            zlims = [min(zs)-0.08, max(zs)+0.08];
        else
            xlims = [-1 1]; ylims = [-1 1]; zlims = [-1 1];
        end

        subplot(1,3,1); hold on
        if any(valid)
            plotTriangulatedCheetah2D(P, varsUsed, [1 2]);
            % scatter(P(valid,1), P(valid,2), 50, 'filled');
            % draw2DSkeleton(P, varsUsed, [1 2]);
            % label2DPoints(P, varsUsed, [1 2]);
        end
        xlabel('X'); ylabel('Y');
        title(sprintf('Top (XY), abs frame %d', absFrames(f)));
        axis equal; grid on; xlim(xlims); ylim(ylims);

        subplot(1,3,2); hold on
        if any(valid)
            plotTriangulatedCheetah2D(P, varsUsed, [1 3]);
            % scatter(P(valid,1), P(valid,3), 50, 'filled');
            % draw2DSkeleton(P, varsUsed, [1 3]);
            % label2DPoints(P, varsUsed, [1 3]);
        end
        xlabel('X'); ylabel('Z');
        title(sprintf('Side (XZ), abs frame %d', absFrames(f)));
        axis equal; grid on; xlim(xlims); ylim(zlims);

        subplot(1,3,3); hold on
        if any(valid)
            plotTriangulatedCheetah2D(P, varsUsed, [2 3]);
            % scatter(P(valid,2), P(valid,3), 50, 'filled');
            % draw2DSkeleton(P, varsUsed, [2 3]);
            % label2DPoints(P, varsUsed, [2 3]);
        end
        xlabel('Y'); ylabel('Z');
        title(sprintf('Front (YZ), abs frame %d', absFrames(f)));
        axis equal; grid on; xlim(ylims); ylim(zlims);

        sgtitle(sprintf('Multiview triangulated orthographic views (%s)', videoTag));
        drawnow
        writeVideo(vw, getframe(fg));
    end

    close(vw);
    close(fg);

    fprintf('Orthographic video saved: %s\n', outName);
end
%[text] 
%[text] ## -------------------- Inspect one frame on raw image --------------------
if doInspect 
    camIdx = find(camIDs == inspectCam, 1);
    assert(~isempty(camIdx), 'inspectCam not in camIDs');

    rawFrameIdx = startI + offsets(camIdx) + (inspectFrame - 1);
    videoPath   = fullfile(parent_folder, sprintf("cam%d_mjpeg.avi", inspectCam));

    if exist(videoPath, 'file') == 2
        vr = VideoReader(videoPath);
        t = (rawFrameIdx - 1) / vr.FrameRate;
        t = max(0, min(t, vr.Duration - 1/vr.FrameRate));
        vr.CurrentTime = t;
        img = readFrame(vr);

        figure('Color','w'); 
        imshow(img); hold on 

        Pframe = squeeze(points3D_raw(inspectFrame,:,:)); % [nBody x 3]
        valid = all(isfinite(Pframe),2);

        % Distorted coords for plotting on the raw video frame
        uv_reproj = nan(nBody,2);
        for b = 1:nBody
            if valid(b)
                uv_reproj(b,:) = world2img(Pframe(b,:), invert(posesPlot{camIdx}), intrinsicParams{camIdx}, ...
                    'ApplyDistortion', true);
            end
        end

        camname = sprintf('cam%d', inspectCam);
        rawTbl  = rawLabelsCropped.(camname);

        uv_manual = nan(nBody,2);
        for b = 1:nBody
            bp = varsUsed{b};
            uv_manual(b,:) = [rawTbl.(sprintf('%s_x',bp))(inspectFrame), ...
                              rawTbl.(sprintf('%s_y',bp))(inspectFrame)];
        end

        plot(uv_manual(:,1), uv_manual(:,2), 'r+', 'MarkerSize', 10, 'LineWidth', 2); 
        plot(uv_reproj(:,1), uv_reproj(:,2), 'go', 'MarkerSize', 8, 'LineWidth', 2); 

        for b = 1:nBody
            if all(isfinite(uv_manual(b,:)))
                text(uv_manual(b,1), uv_manual(b,2), ['  manual-' varsUsed{b}], ... 
                    'Color','r','FontSize',10,'FontWeight','bold'); 
            end
            if all(isfinite(uv_reproj(b,:)))
                text(uv_reproj(b,1), uv_reproj(b,2), ['  reproj-' varsUsed{b}], ...
                    'Color','g','FontSize',10,'FontWeight','bold');
            end
        end

        legend({'manual raw labels','reprojected triangulated 3D'}, 'Location','best'); 
        title(sprintf('Cam%d raw frame %d | multiview triangulated frame %d (abs frame %d)', ... 
            inspectCam, rawFrameIdx, inspectFrame, absFrames(inspectFrame))); 
        hold off 
    else
        warning('Video not found: %s', videoPath);
    end
end 
%%
%[text] ## -------------------- Reprojection overlay: best and worst frame --------------------
%[text] reprojRMSE_cam_perframe (per-frame per-camera in-sample RMSE) is computed
%[text] right after triangulation above and reused here.

overlayVideos = arrayfun(@(c) fullfile(parent_folder, sprintf('cam%d_mjpeg.avi', c)), ...
    camIDs, 'UniformOutput', false);
haveOverlayVideos = doFigs && all(cellfun(@(p) exist(p,'file')==2, overlayVideos));
if doFigs && ~haveOverlayVideos
    warning('multiviewRecon:noVideoForOverlay', ...
        'Skipping reprojection overlays: cam<N>_mjpeg.avi missing in %s', char(parent_folder));
end
if haveOverlayVideos
[~, bestFrame]  = min(reprojRMSE_frame,  [], 'omitnan');  bestFrame  = bestFrame(1);
[~, worstFrame] = max(reprojRMSE_frame,  [], 'omitnan');  worstFrame = worstFrame(1);

fprintf('\nReprojection overlays: best LSQ %d (abs %d, RMSE %.2f px) | worst LSQ %d (abs %d, RMSE %.2f px)\n', ...  
    bestFrame,  startI + bestFrame  - 1, reprojRMSE_frame(bestFrame), ... 
    worstFrame, startI + worstFrame - 1, reprojRMSE_frame(worstFrame));  

partColors = [
    0.0  0.4  1.0;   % l_eye
    1.0  0.2  0.2;   % r_eye
    0.0  0.0  0.0;   % nose
    0.0  0.8  0.8;   % l_ear
    0.8  0.0  0.8;   % r_ear
    0.9  0.7  0.0;   % spine
    0.0  0.7  0.0;   % neck_base
];

ncols_fig = ceil(ncams / 2);
nrows_fig = ceil(ncams / ncols_fig);

frameSet    = [bestFrame,  worstFrame];
frameLabels = {'best',     'worst'};

for fi = 1:numel(frameSet) 
    fIdx   = frameSet(fi);
    fLabel = frameLabels{fi};

    Pframe = squeeze(points3D_video(fIdx,:,:));  % [nBody x 3]

    figure('Color','w', ...
        'Name', sprintf('Reprojection overlay (%s), LSQ frame %d', fLabel, fIdx), ...
        'Visible', onoff(showFigures), ...
        'Position', [60 60 420*ncols_fig 400*nrows_fig]);

    sgtitle(sprintf('Reprojection overlay — %s frame (LSQ idx %d, abs %d, RMSE %.2f px)', ...  
        fLabel, fIdx, startI + fIdx - 1, reprojRMSE_frame(fIdx)), 'FontWeight','bold');  

    for k = 1:ncams
        camNum  = camIDs(k);
        camname = sprintf('cam%d', camNum);

        videoPath = fullfile(parent_folder, sprintf('cam%d_mjpeg.avi', camNum));
        assert(exist(videoPath,'file')==2, 'Video not found: %s', videoPath);

        rawFrameIdx = startI + offsets(k) + (fIdx - 1);
        vr = VideoReader(videoPath);
        t  = max(0, min((rawFrameIdx-1)/vr.FrameRate, vr.Duration - 1/vr.FrameRate));
        vr.CurrentTime = t;
        img = readFrame(vr);

        subplot(nrows_fig, ncols_fig, k);
        imshow(img); hold on;

        for b = 1:nBody
            bp  = varsUsed{b};
            col = partColors(b,:);

            xName = sprintf('%s_x', bp);
            yName = sprintf('%s_y', bp);

            % Undistorted coords — validity gating only
            if all(isfinite(Pframe(b,:)))
                pixels_proj = world2img(Pframe(b,:), invert(posesPlot{k}), intrinsicParams{k});
            else
                pixels_proj = [NaN NaN];
            end
            pixels_obs_ud = [undistortedLabels.(camname).(xName)(fIdx), ...
                             undistortedLabels.(camname).(yName)(fIdx)];
            llh_val = likelihoods.(camname).(bp)(fIdx);

            valid = llh_val > visibilityThreshold && ...
                    all(isfinite(pixels_proj)) && all(isfinite(pixels_obs_ud)) && ...
                    all(pixels_proj > 1) && all(pixels_obs_ud > 1);

            % Distorted coords — plotting on raw video frame
            if all(isfinite(Pframe(b,:)))
                pixels_proj_dist = world2img(Pframe(b,:), invert(posesPlot{k}), intrinsicParams{k}, ...
                    'ApplyDistortion', true);
            else
                pixels_proj_dist = [NaN NaN];
            end
            rx = rawLabelsCropped.(camname).(xName)(fIdx);
            ry = rawLabelsCropped.(camname).(yName)(fIdx);
            if isfinite(rx) && isfinite(ry)
                pixels_obs_raw = [rx ry];
            else
                pixels_obs_raw = [NaN NaN];
            end

            % Always draw projection if 3D point exists
            if all(isfinite(pixels_proj_dist)) && all(pixels_proj_dist > 1)
                plot(pixels_proj_dist(1), pixels_proj_dist(2), 'o', ...
                    'Color', col, 'MarkerSize', 8, 'LineWidth', 2);
                text(pixels_proj_dist(1)+5, pixels_proj_dist(2), bp, ...
                    'Color', col, 'FontSize', 7, 'FontWeight','bold');
            end

            % Draw observation and error line only if LLH-gated valid
            if valid && all(isfinite(pixels_obs_raw))
                plot(pixels_obs_raw(1), pixels_obs_raw(2), '+', ...
                    'Color', col, 'MarkerSize', 10, 'LineWidth', 2);
                plot([pixels_proj_dist(1) pixels_obs_raw(1)], ...
                     [pixels_proj_dist(2) pixels_obs_raw(2)], ...
                     '-', 'Color', col, 'LineWidth', 1.5);
            end
        end

        rmse_cam = reprojRMSE_cam_perframe(fIdx, k);
        title(sprintf('cam%d | raw frame %d | RMSE %.2f px', camNum, rawFrameIdx, rmse_cam), ...
            'Interpreter','none');
        hold off;
    end

    if saveFigures
        savepath = fullfile(reconFolder, ...
            sprintf('Reprojection_overlay_%s_frame%d.png', fLabel, fIdx));
        exportgraphics(gcf, savepath, 'Resolution', 300);
    end
end
end

fprintf('\nDone.\n');
end

% ========================================================================
% Local functions
% ========================================================================
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

function S = packPointsStruct(points3D, varsUsed)
    nBody = numel(varsUsed);
    S = struct();
    for b = 1:nBody
        S.(varsUsed{b}) = squeeze(points3D(:,b,:)); % [nframes x 3]
    end
end

function pointsOut = filter3DTrajectories(pointsIn, modeName, medianWindow, lpfCfg)
    pointsOut = pointsIn;
    [nframes, nBody, ~] = size(pointsIn);

    for b = 1:nBody
        for d = 1:3
            x = squeeze(pointsIn(:,b,d));

            if all(~isfinite(x))
                continue
            end

            valid = isfinite(x);
            t = (1:nframes).';

            xFill = x;
            if sum(valid) >= 2
                xFill(~valid) = interp1(t(valid), x(valid), t(~valid), 'linear', 'extrap');
            else
                continue
            end

            switch lower(modeName)
                case 'median'
                    y = medfilt1(xFill, medianWindow, 'truncate');

                case 'lpf'
                    fc = lpfCfg.fc;
                    fs = lpfCfg.fs;
                    [bb, aa] = butter(4, fc/(fs/2), 'low');
                    y = filtfilt(bb, aa, xFill);

                otherwise
                    error('Unknown filter mode: %s', modeName);
            end

            y(~valid) = nan;
            pointsOut(:,b,d) = y;
        end
    end
end

function plotPointTrajectories(points3D, varsUsed, absFrames, figTitleText)
    [~, nBody, ~] = size(points3D);

    figure('Color','w');
    set(gcf,'Position',[80 80 1600 900]);

    for b = 1:nBody
        subplot(ceil(nBody/2), 2, b);
        X = squeeze(points3D(:,b,1));
        Y = squeeze(points3D(:,b,2));
        Z = squeeze(points3D(:,b,3));

        plot(absFrames, X, 'LineWidth', 1.0); hold on
        plot(absFrames, Y, 'LineWidth', 1.0);
        plot(absFrames, Z, 'LineWidth', 1.0);
        grid on
        xlabel('Absolute frame')
        ylabel('Position (m)')
        title(varsUsed{b}, 'Interpreter','none')
        legend('X','Y','Z', 'Location','best')
    end

    sgtitle(figTitleText);
end

function drawSegmentIfValid(P, varsUsed, a, b, styleStr, lw)
    ia = find(strcmp(varsUsed, a), 1);
    ib = find(strcmp(varsUsed, b), 1);
    if isempty(ia) || isempty(ib), return; end
    if all(isfinite(P(ia,:))) && all(isfinite(P(ib,:)))
        plot3([P(ia,1) P(ib,1)], [P(ia,2) P(ib,2)], [P(ia,3) P(ib,3)], ...
            styleStr, 'LineWidth', lw);
    end
end

function drawMidToPointIfValid(P, varsUsed, a, b, c, styleStr, lw)
    ia = find(strcmp(varsUsed, a), 1);
    ib = find(strcmp(varsUsed, b), 1);
    ic = find(strcmp(varsUsed, c), 1);
    if isempty(ia) || isempty(ib) || isempty(ic), return; end

    if all(isfinite(P(ia,:))) && all(isfinite(P(ib,:))) && all(isfinite(P(ic,:)))
        mid = (P(ia,:) + P(ib,:))/2;
        plot3([mid(1) P(ic,1)], [mid(2) P(ic,2)], [mid(3) P(ic,3)], ...
            styleStr, 'LineWidth', lw);
    end
end

function draw2DSkeleton(P, varsUsed, dims)
    draw2DSegment(P, varsUsed, 'l_eye', 'r_eye', dims, 'k-', 2);
    draw2DMidToPoint(P, varsUsed, 'l_eye', 'r_eye', 'nose', dims, 'k-', 2);
    draw2DSegment(P, varsUsed, 'neck_base', 'spine', dims, 'k-', 2);
    draw2DSegment(P, varsUsed, 'l_ear', 'r_ear', dims, 'k-', 2);
    draw2DSegment(P, varsUsed, 'l_ear', 'l_eye', dims, 'k-', 2);
    draw2DSegment(P, varsUsed, 'r_ear', 'r_eye', dims, 'k-', 2);
end

function draw2DSegment(P, varsUsed, a, b, dims, styleStr, lw)
    ia = find(strcmp(varsUsed, a), 1);
    ib = find(strcmp(varsUsed, b), 1);
    if isempty(ia) || isempty(ib), return; end

    if all(isfinite(P(ia,:))) && all(isfinite(P(ib,:)))
        plot([P(ia,dims(1)) P(ib,dims(1))], [P(ia,dims(2)) P(ib,dims(2))], ...
            styleStr, 'LineWidth', lw);
    end
end

function draw2DMidToPoint(P, varsUsed, a, b, c, dims, styleStr, lw)
    ia = find(strcmp(varsUsed, a), 1);
    ib = find(strcmp(varsUsed, b), 1);
    ic = find(strcmp(varsUsed, c), 1);
    if isempty(ia) || isempty(ib) || isempty(ic), return; end

    if all(isfinite(P(ia,:))) && all(isfinite(P(ib,:))) && all(isfinite(P(ic,:)))
        mid = (P(ia,:) + P(ib,:))/2;
        plot([mid(dims(1)) P(ic,dims(1))], [mid(dims(2)) P(ic,dims(2))], ...
            styleStr, 'LineWidth', lw);
    end
end

function label2DPoints(P, varsUsed, dims)
    for i = 1:size(P,1)
        if all(isfinite(P(i,:)))
            text(P(i,dims(1)), P(i,dims(2)), ['  ' varsUsed{i}], ...
                'FontSize', 9, 'Color', 'k');
        end
    end
end
