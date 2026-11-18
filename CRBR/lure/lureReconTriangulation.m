function lureReconTriangulation(parent_folder, camIDs, opts)
% lureReconTriangulation  CRBR lure triangulation from AcinoSet image annotations.
%
%   lureReconTriangulation(parent_folder, camIDs, opts)
%
% Direct multi-view triangulation (triangulateMultiview) of the lure marker
% from Blender/manual labels. Follows the same camera,
% offset, frame-range and undistortion conventions as the reconstruction
% pipeline. Writes the lure 3D track + reprojection diagnostics to
% <parent_folder>/lureRecon/.
%
% opts fields (all optional; defaults are headless + fast for batch):
%   lureName('lure')  visibilityThreshold(0.5)  minViewsToTriangulate(2)
%   reprojOutlierThreshPx(10)  plotCameraConfig(false)
%   lureFrameRange([]) useCheetahLimits(true)
%   outputFileName('3DLurePos.mat') savePerCameraLabels(true)
%   showFigures(false)  saveFigures(true)
%
% Figure show/save convention (shared across the pipeline):
%   showFigures - pop up the figure windows
%   saveFigures - write the report PNGs into <parent_folder>/lureRecon/
% By default the reprojection-quality summary (lure_reprojection_quality.png)
% and the best-RMSE overlay (lure_reproj_best_frame.png) are SAVED but NOT
% shown: the reprojection check is preserved for auditing without interrupting
% batch runs.
%
% Called with no arguments it runs a default dataset (standalone debugging).

if nargin < 3 || isempty(opts)
    opts = struct();
end
if nargin < 2 || isempty(camIDs)
    camIDs = [1 2 3 4 5 6];
end
if nargin < 1 || isempty(parent_folder)
    datasetName   = "phantomRun_20190307";
    repoRoot      = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    parent_folder = resolveDatasetFolder(datasetName, repoRoot);
else
    [~, nm]     = fileparts(char(parent_folder));
    datasetName = string(nm);
end

ncameras = numel(camIDs);

lureName                  = getOpt(opts, 'lureName', 'lure');
visibilityThreshold       = getOpt(opts, 'visibilityThreshold', 0.5);   % same gate as cheetahError
minViewsToTriangulate     = getOpt(opts, 'minViewsToTriangulate', 2);
REPROJ_OUTLIER_THRESH_PX  = getOpt(opts, 'reprojOutlierThreshPx', 10);
frameRange                = getOpt(opts, 'lureFrameRange', []);
useCheetahLimits          = getOpt(opts, 'useCheetahLimits', true);
outputFileName            = char(getOpt(opts, 'outputFileName', '3DLurePos.mat'));
savePerCameraLabels       = getOpt(opts, 'savePerCameraLabels', true);

sf                        = getOpt(opts, 'showFigures', false);
saveFigures               = getOpt(opts, 'saveFigures', true);
PLOT_CAMERA_CONFIG        = getOpt(opts, 'plotCameraConfig',    sf);
doReprojFigs              = sf || saveFigures;

outDir = fullfile(parent_folder, "lureRecon");

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf("Dataset: %s\n", datasetName);  
fprintf("Folder : %s\n", parent_folder);  
fprintf("Cameras: %s\n", mat2str(camIDs));  
%%
%[text] \-------------------- Load labels + LLH --------------------
labelTables = cell(1, ncameras);
llhTables   = cell(1, ncameras);

for ii = 1:ncameras  
    cam = camIDs(ii);
    labelPath = fullfile(parent_folder, sprintf("cam%d_labels_wide.csv", cam));
    llhPath   = fullfile(parent_folder, sprintf("cam%d_llh_wide.csv", cam));

    assert(exist(labelPath, 'file') == 2, "Missing label file: %s", labelPath);
    assert(exist(llhPath, 'file') == 2, "Missing LLH file: %s", llhPath);

    labelTables{ii} = readtable(labelPath, VariableNamingRule="preserve");
    llhTables{ii}   = readtable(llhPath, VariableNamingRule="preserve");

    fprintf("cam%d: labels %s, llh %s\n", ...  
        cam, mat2str(size(labelTables{ii})), mat2str(size(llhTables{ii})));  
end  
%%
%[text] \-------------------- Offsets --------------------
offsetPath = dir(fullfile(parent_folder, "AcinoSet_data", "*offset*"));
if ~isempty(offsetPath)  
    offsets = readtable(fullfile(offsetPath(1).folder, offsetPath(1).name), ...
        VariableNamingRule="preserve");
    fprintf("Loaded offsets from %s\n", offsetPath(1).name);  
else
    offsets = table(camIDs(:), zeros(ncameras, 1), ...
        'VariableNames', {'cam', 'offset'});
    fprintf("No offset file found. Using zero offsets.\n");
end  

assert(ismember("offset", string(offsets.Properties.VariableNames)), ...
    "Offsets table must contain an offset column.");

for ii = 1:ncameras  
    fprintf("  cam%d offset = %g\n", camIDs(ii), getCamOffset(offsets, camIDs(ii)));  
end  
%%
%[text] \-------------------- Load camera parameters --------------------
jsonPath = resolveCamSceneSba(parent_folder, ncameras);
fprintf("Calibration: %s\n", jsonPath);  

extrinsics = jsondecode(fileread(jsonPath));
imageSize  = [extrinsics.camera_resolution(2), extrinsics.camera_resolution(1)];

poses = cell(1, ncameras);
intrinsicParams = cell(1, ncameras);
cameraPoses = table( ...
    uint32(camIDs(:)), ...
    repmat(rigidtform3d, ncameras, 1), ...
    'VariableNames', {'ViewId', 'AbsolutePose'});

for ii = 1:ncameras
    c = extrinsics.cameras(ii);

    R_cw = c.r;       % world -> camera
    t_cw = c.t(:);    % 3x1

    R_wc = R_cw.';
    C_w  = -(R_wc * t_cw);

    T_wc = rigidtform3d(R_wc, C_w.');
    poses{ii} = T_wc;
    cameraPoses.AbsolutePose(ii) = T_wc;

    intrinsicParams{ii} = cameraIntrinsicsFromOpenCV(c.k, c.d, imageSize);
end

intrinsicsForTriangulation = [intrinsicParams{:}];

if PLOT_CAMERA_CONFIG  
    figure('Color', 'w', 'Name', 'Camera configuration');  
    if exist("plotCameraConfig6", "file") == 2 && ncameras == 6
        plotCameraConfig6(poses, 0.1);  
    elseif exist("plotCameraConfigLT6", "file") == 2
        plotCameraConfigLT6(poses, 0.1, camIDs);
    end
    title("Camera configuration");  
    axis vis3d  
    grid on  
end  
%%
%[text] \-------------------- Reconstruction frame range --------------------
refLabels = labelTables{1};
limitFile = dir(fullfile(parent_folder, "*limit*"));

if ~isempty(frameRange)
    assert(numel(frameRange)==2 && all(isfinite(frameRange)), ...
        'lureFrameRange must be [start end].');
    startI=round(frameRange(1)); endI=round(frameRange(2));
    fprintf("Using explicit lure-only reconstruction range: absolute frame %d -> %d\n",startI,endI);
elseif useCheetahLimits && ~isempty(limitFile)
    lims = readtable(fullfile(limitFile(1).folder, limitFile(1).name), ...
        VariableNamingRule="preserve");
    startI = lims{1, 1};
    endI   = lims{1, 2};
    fprintf("Using reconstruction range from limit file: absolute frame %d -> %d\n", ...  
        startI, endI);  
else
    offVec=arrayfun(@(cam)getCamOffset(offsets,cam),camIDs);
    heights=cellfun(@height,labelTables);
    startI=max([1,1-offVec]);
    endI=min(heights-offVec);
    fprintf("Using full lure-label range independent of cheetah limits: absolute frame %d -> %d\n", ...
        startI, endI);
end  

startI = max(1, startI);
offVec=arrayfun(@(cam)getCamOffset(offsets,cam),camIDs);
heights=cellfun(@height,labelTables);
endI=min([endI,heights-offVec]);
assert(startI <= endI, "Invalid reconstruction frame range.");

absFrames = (startI:endI).';
nframes   = numel(absFrames);
fprintf("Final reconstruction range: %d -> %d (nframes=%d)\n", ...    
    startI, endI, nframes);    
%%
%[text] \-------------------- Crop + undistort lure labels --------------------
rawLabels    = struct();
undistorted  = struct();
likelihoods  = struct();

xCol = [lureName '_x'];
yCol = [lureName '_y'];
lCol = [lureName '_LLH'];

for ii = 1:ncameras  
    cam = camIDs(ii);
    camname = sprintf('cam%d', cam);
    offset = getCamOffset(offsets, cam);

    labelsTbl = labelTables{ii};
    llhTbl    = llhTables{ii};

    assert(ismember(xCol, labelsTbl.Properties.VariableNames), ...
        "cam%d is missing %s", cam, xCol);
    assert(ismember(yCol, labelsTbl.Properties.VariableNames), ...
        "cam%d is missing %s", cam, yCol);
    assert(ismember(lCol, llhTbl.Properties.VariableNames), ...
        "cam%d is missing %s", cam, lCol);

    rowStart = startI + offset;
    rowEnd   = endI + offset;

    assert(rowStart >= 1, "cam%d: rowStart=%d < 1", cam, rowStart);
    assert(rowEnd <= height(labelsTbl), ...
        "cam%d: rowEnd=%d exceeds label rows=%d", cam, rowEnd, height(labelsTbl));
    assert(rowEnd <= height(llhTbl), ...
        "cam%d: rowEnd=%d exceeds LLH rows=%d", cam, rowEnd, height(llhTbl));

    labelsCrop = labelsTbl(rowStart:rowEnd, :);
    llhCrop    = llhTbl(rowStart:rowEnd, :);

    rawXY = [labelsCrop.(xCol), labelsCrop.(yCol)];
    llh   = llhCrop.(lCol);

    validRaw = all(isfinite(rawXY), 2) & isfinite(llh) & ...
        llh > visibilityThreshold & all(rawXY > 1, 2);

    xyUD = nan(nframes, 2);
    if any(validRaw)
        xyUD(validRaw, :) = undistortPoints(rawXY(validRaw, :), intrinsicParams{ii});
    end

    rawLabels.(camname).XY       = rawXY;
    rawLabels.(camname).Index    = absFrames;
    rawLabels.(camname).LLH      = llh;
    rawLabels.(camname).Valid    = validRaw;
    rawLabels.(camname).RawFrame = absFrames + offset;

    undistorted.(camname).XY    = xyUD;
    undistorted.(camname).Index = absFrames;
    likelihoods.(camname)       = llh;

    results = table(absFrames(validRaw), rawXY(validRaw, 1), rawXY(validRaw, 2), ...
        'VariableNames', {'Index', 'X', 'Y'});
    if savePerCameraLabels
        save(fullfile(parent_folder, sprintf("cam%dlureLabels.mat", cam)), "results");
    end

    fprintf("cam%d: %d/%d lure observations after LLH/finite gate", ...            
        cam, height(results), nframes);            
    if height(results) > 0
        fprintf(" (abs %d -> %d)", min(results.Index), max(results.Index));            
    end
    fprintf("\n");

end  
%%
%[text] \-------------------- Multi-view triangulation --------------------
fprintf("\n=== Multi-view triangulating lure ===\n");  

Lure3D = struct();
Lure3D.Index          = [];
Lure3D.Mark           = {};
Lure3D.IsInterpolated = [];
Lure3D.NCams          = [];
Lure3D.ViewIds        = {};

triangulateReturnedErr = [];
validIndex             = [];

kOut = 1;
for fi = 1:nframes
    viewIDs = zeros(0, 1, 'uint32');
    pts2D   = zeros(0, 2);

    for ii = 1:ncameras
        cam = camIDs(ii);
        camname = sprintf('cam%d', cam);

        uv  = undistorted.(camname).XY(fi, :);
        llh = likelihoods.(camname)(fi);

        if all(isfinite(uv)) && isfinite(llh) && llh > visibilityThreshold
            viewIDs(end+1, 1) = uint32(cam); 
            pts2D(end+1, :)   = uv;          
        end
    end

    if numel(viewIDs) < minViewsToTriangulate
        continue
    end

    track = pointTrack(viewIDs.', pts2D);

    try
        [X, reprojErr, validIdx] = triangulateMultiview( ...
            track, cameraPoses, intrinsicsForTriangulation);
    catch ME
        warning("triangulateMultiview failed at absolute frame %d: %s", ...
            absFrames(fi), ME.message);
        continue
    end

    if isempty(X) || ~all(isfinite(X(1, :))) || isempty(validIdx) || ~validIdx(1)
        continue
    end

    Lure3D.Index(kOut)          = absFrames(fi);
    Lure3D.Mark{kOut}           = X(1, :);
    Lure3D.IsInterpolated(kOut) = false;
    Lure3D.NCams(kOut)          = numel(viewIDs);
    Lure3D.ViewIds{kOut}        = viewIDs.';

    triangulateReturnedErr(kOut, 1) = reprojErr(1); 
    validIndex(kOut, 1)             = validIdx(1);  
    kOut = kOut + 1;
end

Lure3D.TriangulateReprojErr_px = triangulateReturnedErr.';
Lure3D.ValidIndex              = validIndex.';

fprintf("Triangulated %d real 3D lure points.\n", numel(Lure3D.Index));  
printCoverageSummary(Lure3D, undistorted, camIDs);  
%%
%[text] \-------------------- Reprojection validation -------------------- CRBR convention: RMSE is computed in undistorted image space, while the visual overlay projects distorted pixels onto the raw video frame.
[reprojErr_frame, reprojErr_camPx] = computeReprojectionErrors( ...
    Lure3D, undistorted, poses, intrinsicParams, camIDs, startI);

Lure3D.ReprojErr_px    = reprojErr_frame.';
Lure3D.ReprojErr_camPx = reprojErr_camPx;

printReprojectionSummary( ...    
    Lure3D.Index, reprojErr_frame, reprojErr_camPx, camIDs, REPROJ_OUTLIER_THRESH_PX);    

% Reprojection-quality summary (lure_reprojection_quality.png) and the
% best-RMSE overlay (lure_reproj_best_frame.png). Saved when saveFigures is on
% (default), displayed only when showFigures is on.
if doReprojFigs
    plotReprojectionQuality( ...
        Lure3D, reprojErr_frame, reprojErr_camPx, camIDs, ...
        REPROJ_OUTLIER_THRESH_PX, outDir, sf, saveFigures);

    writeBestFrameOverlay( ...
        Lure3D, reprojErr_frame, reprojErr_camPx, labelTables, llhTables, ...
        offsets, startI, poses, intrinsicParams, camIDs, parent_folder, ...
        outDir, lureName, visibilityThreshold, sf, saveFigures);
end

save(fullfile(parent_folder, outputFileName), "Lure3D");
fprintf("Saved %s\n", fullfile(parent_folder, outputFileName));
end

% ======================================================================== Helpers ========================================================================
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

function parentFolder = resolveDatasetFolder(datasetName, repoRoot)
    candidates = [
        fullfile(repoRoot, "outputs", datasetName)
        fullfile(pwd, "outputs", datasetName)
        fullfile(fileparts(pwd), "outputs", datasetName)
        fullfile("outputs", datasetName)
    ];

    for i = 1:numel(candidates)
        if isfolder(candidates(i))
            parentFolder = candidates(i);
            return
        end
    end

    error("Could not find data folder for dataset %s.", datasetName);
end

function jsonPath = resolveCamSceneSba(parentFolder, ncameras)
    acinoDir = fullfile(parentFolder, "AcinoSet_data");

    exactPath = fullfile(acinoDir, sprintf("%d_cam_scene_sba.json", ncameras));
    if isfile(exactPath)
        jsonPath = exactPath;
        return
    end

    sharedPath = fullfile(acinoDir, "cam_scene_sba.json");
    if isfile(sharedPath)
        jsonPath = sharedPath;
        return
    end

    matches = dir(fullfile(acinoDir, "*cam_scene_sba.json"));
    assert(~isempty(matches), "No *cam_scene_sba.json found in %s.", acinoDir);
    jsonPath = fullfile(matches(1).folder, matches(1).name);
end

function offset = getCamOffset(offsets, cam)
    names = string(offsets.Properties.VariableNames);

    if ismember("cam", names)
        idCol = "cam";
    elseif ismember("camID", names)
        idCol = "camID";
    else
        idCol = "";
    end

    if idCol ~= ""
        row = find(offsets.(idCol) == cam, 1);
        assert(~isempty(row), "Camera %d not found in offsets table.", cam);
        offset = offsets.offset(row);
    else
        assert(height(offsets) >= cam, ...
            "Offsets table has no cam/camID column and fewer than cam%d rows.", cam);
        offset = offsets.offset(cam);
    end

    if isempty(offset) || ~isfinite(offset)
        offset = 0;
    end
    offset = double(offset);
end

function printCoverageSummary(Lure3D, undistorted, camIDs)
    if isempty(Lure3D.Index)
        fprintf("No frames had enough camera observations for triangulation.\n");
        return
    end

    firstTriAbs = min(Lure3D.Index);
    lastTriAbs  = max(Lure3D.Index);
    nPossible   = lastTriAbs - firstTriAbs + 1;
    nMissing    = nPossible - numel(Lure3D.Index);

    fprintf("  Absolute frame range : %d -> %d\n", firstTriAbs, lastTriAbs);
    fprintf("  Coverage             : %d / %d frames (%.1f%%)\n", ...
        numel(Lure3D.Index), nPossible, 100 * numel(Lure3D.Index) / nPossible);
    fprintf("  Missing frames       : %d\n", nMissing);
    fprintf("  Cameras per frame    : avg %.1f, min %d, max %d\n", ...
        mean(Lure3D.NCams), min(Lure3D.NCams), max(Lure3D.NCams));

    for ii = 1:numel(camIDs)
        cam = camIDs(ii);
        camname = sprintf('cam%d', cam);
        valid = all(isfinite(undistorted.(camname).XY), 2);
        nContrib = sum(ismember(Lure3D.Index, undistorted.(camname).Index(valid)));
        fprintf("    cam%d: visible in %d / %d triangulated frames\n", ...
            cam, nContrib, numel(Lure3D.Index));
    end
end

function [reprojErr_frame, reprojErr_camPx] = computeReprojectionErrors( ...
        Lure3D, undistorted, poses, intrinsicParams, camIDs, startI)

    nTri = numel(Lure3D.Index);
    ncameras = numel(camIDs);
    reprojErr_camPx = nan(nTri, ncameras);

    for ki = 1:nTri
        absFrame = Lure3D.Index(ki);
        fidx = absFrame - startI + 1;
        if fidx < 1
            continue
        end

        pt3D = Lure3D.Mark{ki}(:).';
        usedViews = Lure3D.ViewIds{ki};

        for ii = 1:ncameras
            cam = camIDs(ii);
            if ~ismember(uint32(cam), uint32(usedViews))
                continue
            end

            camname = sprintf('cam%d', cam);
            if fidx > size(undistorted.(camname).XY, 1)
                continue
            end

            obs2D = undistorted.(camname).XY(fidx, :);
            proj2D = world2img(pt3D, invert(poses{ii}), intrinsicParams{ii});

            if all(isfinite(obs2D)) && all(isfinite(proj2D)) && ...
                    all(obs2D > 1) && all(proj2D > 1)
                reprojErr_camPx(ki, ii) = norm(proj2D - obs2D);
            end
        end
    end

    reprojErr_frame = nan(nTri, 1);
    for ki = 1:nTri
        vals = reprojErr_camPx(ki, isfinite(reprojErr_camPx(ki, :)));
        if ~isempty(vals)
            reprojErr_frame(ki) = sqrt(mean(vals.^2));
        end
    end
end

function printReprojectionSummary(absIdx, reprojErr_frame, reprojErr_camPx, camIDs, outlierThresh)
    fprintf("\n=== Reprojection validation (saved triangulated lure) ===\n");

    validF = isfinite(reprojErr_frame);
    if ~any(validF)
        fprintf("  No valid reprojection errors to summarize.\n");
        return
    end

    fprintf("  Per-frame RMSE (px): mean=%.2f  median=%.2f  p95=%.2f  max=%.2f\n", ...
        mean(reprojErr_frame(validF)), median(reprojErr_frame(validF)), ...
        prctile(reprojErr_frame(validF), 95), max(reprojErr_frame(validF)));

    fprintf("  Per-camera RMSE (px):\n");
    for ii = 1:numel(camIDs)
        vals = reprojErr_camPx(:, ii);
        vals = vals(isfinite(vals));
        if isempty(vals)
            fprintf("    cam%d: no observations\n", camIDs(ii));
        else
            fprintf("    cam%d: mean=%.2f  median=%.2f  p95=%.2f  max=%.2f\n", ...
                camIDs(ii), sqrt(mean(vals.^2)), median(vals), prctile(vals, 95), max(vals));
        end
    end

    outlierMask = validF & reprojErr_frame > outlierThresh;
    if any(outlierMask)
        fprintf("  WARNING: %d frame(s) with RMSE > %g px:\n", ...
            sum(outlierMask), outlierThresh);
        oidx = find(outlierMask);
        for oi = 1:min(numel(oidx), 10)
            fprintf("    abs frame %d: %.2f px\n", absIdx(oidx(oi)), reprojErr_frame(oidx(oi)));
        end
        if numel(oidx) > 10
            fprintf("    ... and %d more\n", numel(oidx) - 10);
        end
    else
        fprintf("  No outlier frames above %g px threshold.\n", outlierThresh);
    end
end

function plotReprojectionQuality(Lure3D, reprojErr_frame, reprojErr_camPx, camIDs, outlierThresh, outDir, showFig, saveFig)
    validF = isfinite(reprojErr_frame);
    if ~any(validF)
        return
    end

    ncameras = numel(camIDs);
    absIdx = Lure3D.Index(:);
    camColors = lines(ncameras);

    fig = figure('Color', 'w', 'Name', 'Lure reprojection quality', ...
        'Visible', onoff(showFig), 'Position', [100 100 1200 700]);

    subplot(2, 2, 1);
    plot(absIdx, reprojErr_frame, 'k.-', 'MarkerSize', 8, 'LineWidth', 1.2);
    hold on;
    yline(mean(reprojErr_frame(validF)), 'r--', 'Mean', 'LineWidth', 1.2);
    yline(outlierThresh, 'b:', sprintf('%g px', outlierThresh), 'LineWidth', 1);
    grid on; xlabel('Absolute frame'); ylabel('RMSE (px)');
    title('Per-frame reprojection RMSE');

    subplot(2, 2, 2);
    hold on;
    for ii = 1:ncameras
        col = reprojErr_camPx(:, ii);
        validCam = isfinite(col);
        if any(validCam)
            plot(absIdx(validCam), col(validCam), '.-', ...
                'Color', camColors(ii, :), 'MarkerSize', 6, 'LineWidth', 1.0);
        end
    end
    grid on; xlabel('Absolute frame'); ylabel('Error (px)');
    title('Per-camera reprojection error');
    legend(arrayfun(@(c) sprintf("cam%d", c), camIDs, "UniformOutput", false), ...
        'Location', 'best');

    subplot(2, 2, 3);
    bar(absIdx, Lure3D.NCams, 'FaceColor', [0.3 0.5 0.8]);
    grid on; xlabel('Absolute frame'); ylabel('# cameras');
    title('Camera coverage per frame');
    ylim([0, ncameras + 0.5]);

    subplot(2, 2, 4);
    histogram(reprojErr_frame(validF), 'BinWidth', 0.5, ...
        'FaceColor', [0.3 0.6 0.3], 'EdgeColor', 'none');
    hold on;
    xline(mean(reprojErr_frame(validF)), 'r--', 'Mean', 'LineWidth', 1.2);
    xline(prctile(reprojErr_frame(validF), 95), 'b--', 'p95', 'LineWidth', 1.2);
    grid on; xlabel('RMSE (px)'); ylabel('Count');
    title('Distribution of per-frame reprojection RMSE');

    sgtitle('Lure reprojection validation - undistorted pixels');

    if saveFig
        outPath = fullfile(outDir, "lure_reprojection_quality.png");
        print(fig, outPath, '-dpng', '-r150');
        fprintf("  Saved: %s\n", outPath);
    end
end

function writeBestFrameOverlay(Lure3D, reprojErr_frame, reprojErr_camPx, ...
        labelTables, llhTables, offsets, startI, poses, intrinsicParams, ...
        camIDs, parentFolder, outDir, lureName, visibilityThreshold, showFig, saveFig)

    valid = isfinite(reprojErr_frame);
    if ~any(valid)
        return
    end

    validIdx = find(valid);
    [~, bestRel] = min(reprojErr_frame(valid));

    % Only the best-RMSE frame is exported (lure_reproj_best_frame.png).
    overlayFrames = struct( ...
        'ki', {validIdx(bestRel)}, ...
        'label', {'best'});

    ncameras = numel(camIDs);
    ncols = ceil(sqrt(ncameras));
    nrows = ceil(ncameras / ncols);
    xName = sprintf('%s_x', lureName);
    yName = sprintf('%s_y', lureName);
    lName = sprintf('%s_LLH', lureName);

    for si = 1:numel(overlayFrames)
        ki = overlayFrames(si).ki;
        label = overlayFrames(si).label;

        absFrame = Lure3D.Index(ki);
        fIdx = absFrame - startI + 1;
        pt3D = Lure3D.Mark{ki}(:).';
        frameRMSE = reprojErr_frame(ki);

        fig = figure('Color', 'w', ...
            'Name', sprintf('Lure reproj %s frame abs=%d', label, absFrame), ...
            'Visible', onoff(showFig), 'Position', [50 50 400*ncols 320*nrows]);

        for ii = 1:ncameras
            cam = camIDs(ii);
            offset = getCamOffset(offsets, cam);
            rawFrameIdx = startI + offset + (fIdx - 1);

            ax = subplot(nrows, ncols, ii);
            hold(ax, 'on');

            img = readVideoFrameSafe(parentFolder, cam, rawFrameIdx);
            imshow(img, 'Parent', ax);

            % Distorted projection is the coordinate system of the raw frame.
            proj2D = world2img(pt3D, invert(poses{ii}), intrinsicParams{ii}, ...
                'ApplyDistortion', true);

            [label2D, labelVisible, csvFrame] = getRawLureObservation( ...
                labelTables{ii}, llhTables{ii}, rawFrameIdx, ...
                xName, yName, lName, visibilityThreshold);

            if all(isfinite(proj2D))
                plot(ax, proj2D(1), proj2D(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
            end
            if labelVisible
                plot(ax, label2D(1), label2D(2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
            end
            if all(isfinite(proj2D)) && labelVisible
                plot(ax, [proj2D(1) label2D(1)], [proj2D(2) label2D(2)], ...
                    'y-', 'LineWidth', 1.5);
            end

            errStr = "";
            if isfinite(reprojErr_camPx(ki, ii))
                errStr = sprintf(" | %.1f px undist", reprojErr_camPx(ki, ii));
            end

            if isfinite(csvFrame)
                frameStr = sprintf("raw row %g csv %g", rawFrameIdx, csvFrame);
            else
                frameStr = sprintf("raw row %g", rawFrameIdx);
            end

            title(ax, sprintf("cam%d | %s%s", cam, frameStr, errStr), ...
                'FontSize', 9);
        end

        sgtitle(sprintf( ...
            'Lure reprojection - %s RMSE frame (abs %d, RMSE %.2f px undist)', ...
            label, absFrame, frameRMSE));

        if saveFig
            outPath = fullfile(outDir, sprintf("lure_reproj_%s_frame.png", label));
            print(fig, outPath, '-dpng', '-r150');
            fprintf("  Saved overlay (%s RMSE frame, abs=%d): %s\n", ...
                label, absFrame, outPath);
        end
    end
end

function [label2D, visible, csvFrame] = getRawLureObservation( ...
        labelsTbl, llhTbl, rawFrameIdx, xName, yName, lName, visibilityThreshold)

    label2D = [NaN NaN];
    visible = false;
    csvFrame = NaN;

    rawRow = round(rawFrameIdx);
    if abs(rawRow - rawFrameIdx) > eps(rawFrameIdx) || ...
            rawRow < 1 || rawRow > height(labelsTbl) || rawRow > height(llhTbl)
        return
    end

    if ismember("frame", string(labelsTbl.Properties.VariableNames))
        csvFrame = labelsTbl.frame(rawRow);
    end

    if ~ismember(xName, labelsTbl.Properties.VariableNames) || ...
            ~ismember(yName, labelsTbl.Properties.VariableNames) || ...
            ~ismember(lName, llhTbl.Properties.VariableNames)
        return
    end

    label2D = [labelsTbl.(xName)(rawRow), labelsTbl.(yName)(rawRow)];
    llh = llhTbl.(lName)(rawRow);
    visible = all(isfinite(label2D)) && isfinite(llh) && ...
        llh > visibilityThreshold && all(label2D > 1);
end

function img = readVideoFrameSafe(parentFolder, cam, frameIdx)
    try
        videoPath = fullfile(parentFolder, sprintf('cam%d_mjpeg.avi', cam));
        if exist(videoPath, 'file') ~= 2
            videoPath = fullfile(parentFolder, sprintf('cam%d.mp4', cam));
        end
        assert(exist(videoPath, 'file') == 2, 'Video not found: %s', videoPath);

        vr = VideoReader(videoPath);
        if frameIdx < 1 || frameIdx > vr.NumFrames
            img = uint8(30 * ones(108, 192, 3));
            return
        end

        vr.CurrentTime = max(0, min((frameIdx - 1) / vr.FrameRate, ...
            vr.Duration - 1 / vr.FrameRate));
        img = readFrame(vr);
    catch
        img = uint8(30 * ones(108, 192, 3));
    end
end
