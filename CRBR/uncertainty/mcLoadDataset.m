function D = mcLoadDataset(parent_folder, camIDs)
% MCLOADDATASET  Load everything the Monte-Carlo reconstruction needs for one
% dataset, replicating the data ingestion in CRBR reconstruction so the MC fit sees
% identical labels, calibration, scales and frame range.
%
%   D = mcLoadDataset(parent_folder, camIDs)
%
% Returns struct D with fields:
%   poses, intrin      1xncams cells (camera->world tform, cameraIntrinsics)
%   labels0            [nframes x nparts x ncams x 2] UNDISTORTED 2D labels
%   llh                [nframes x nparts x ncams] per-point likelihood
%   scaleHead0, scaleTorso0   nominal rigid-body scales
%   startI, endI, absFrames, nframes  reconstruction frame range (absolute Blender)
%   camIDs, partNames, nparts, ncams
%   sigmas             1xnparts per-part pixel sigma (label noise model)
%   xBarStore          [nframes x 12] nominal poses (per-frame fit seeds)

    ncams     = numel(camIDs);
    partNames = {'l_eye','r_eye','nose','l_ear','r_ear','spine','neck_base'};
    varnames  = partNames;                 % undistortBlenderLabels uses <bp>_x/_y
    nparts    = numel(partNames);

    % -------------------- Rigid-body scales --------------------
    scaleFile = fullfile(parent_folder, 'cheetah_scales.mat');
    if isfile(scaleFile)
        S = load(scaleFile, 'scaleHead', 'scaleTorsoLength');
        D.scaleHead0  = S.scaleHead;
        D.scaleTorso0 = S.scaleTorsoLength;
    else
        warning('mcLoadDataset:noScales', ...
            'cheetah_scales.mat not found in %s; using default scaleHead=scaleTorsoLength=1.', ...
            parent_folder);
        D.scaleHead0  = 1;
        D.scaleTorso0 = 1;
    end

    % -------------------- Camera calibration --------------------
    jf = dir(fullfile(parent_folder, 'AcinoSet_data', '*cam_scene_sba.json'));
    wantedName = sprintf('%d_cam_scene_sba.json', ncams);
    m = find(strcmp({jf.name}, wantedName), 1);
    assert(~isempty(m), 'Could not find %s in AcinoSet_data.', wantedName);
    ext = jsondecode(fileread(fullfile(parent_folder, 'AcinoSet_data', jf(m).name)));

    poses  = cell(1, ncams);
    intrin = cell(1, ncams);
    for ii = 1:ncams
        c    = ext.cameras(ii);
        R_cw = c.r;                 % world -> camera
        t_cw = c.t';
        R_wc = R_cw.';
        C_w  = -(R_cw.' * t_cw.').';
        poses{ii} = rigidtform3d(R_wc, C_w);
        imageSize = [ext.camera_resolution(2) ext.camera_resolution(1)];
        intrin{ii} = cameraIntrinsicsFromOpenCV(c.k, c.d, imageSize);
    end
    D.poses  = poses;
    D.intrin = intrin;

    % -------------------- Offsets (per-camera frame sync) --------------------
    O = readtable(fullfile(parent_folder, 'AcinoSet_data', 'offset.xlsx'), ...
                  'VariableNamingRule', 'preserve');

    % -------------------- Reconstruction frame range --------------------
    ref = readtable(fullfile(parent_folder, sprintf('cam%d_labels_wide.csv', camIDs(1))));
    if ismember('frame', ref.Properties.VariableNames)
        ref = removevars(ref, 'frame');
    end

    limitfile = dir(fullfile(parent_folder, '*limit*'));
    if ~isempty(limitfile)
        lims   = readtable(fullfile(parent_folder, limitfile(1).name), ...
                           'VariableNamingRule', 'preserve');
        startI = lims{1,1};
        endI   = lims{1,2};
    else
        startI = 1;
        endI   = height(ref);
    end
    startI  = max(1, startI);
    endI    = min(height(ref), endI);
    assert(startI <= endI, 'Invalid reconstruction range.');
    expectedAbsFrames = (startI:endI).';

    % -------------------- Nominal poses + absolute-frame provenance --------------------
    optimFile = fullfile(parent_folder, 'storedOptimRaw.mat');
    SS = load(optimFile, 'xBarStore', 'absFrames', 'reconMeta');
    assert(isfield(SS, 'xBarStore'), 'xBarStore not found in %s.', optimFile);
    D.xBarStore = SS.xBarStore;
    if isfield(SS, 'absFrames') && ~isempty(SS.absFrames)
        storedAbsFrames = double(SS.absFrames(:));
        assert(numel(storedAbsFrames) == size(D.xBarStore,1), ...
            '%s contains %d pose rows but %d absFrames entries.', ...
            optimFile, size(D.xBarStore,1), numel(storedAbsFrames));
        assert(isequal(storedAbsFrames, expectedAbsFrames), ...
            ['%s covers absolute frames %d..%d (%d rows), but current limits declare ' ...
             '%d..%d (%d rows). Re-run reconstruction before Monte Carlo.'], ...
            optimFile, storedAbsFrames(1), storedAbsFrames(end), numel(storedAbsFrames), ...
            expectedAbsFrames(1), expectedAbsFrames(end), numel(expectedAbsFrames));
    else
        assert(size(D.xBarStore,1) == numel(expectedAbsFrames), ...
            ['Legacy %s has %d pose rows but current limits declare %d frames (%d..%d). ' ...
             'Re-run reconstruction so absFrames provenance is saved.'], ...
            optimFile, size(D.xBarStore,1), numel(expectedAbsFrames), ...
            expectedAbsFrames(1), expectedAbsFrames(end));
        warning('mcLoadDataset:legacyAbsFrames', ...
            '%s has no absFrames provenance; inferred frames from limits.xlsx.', optimFile);
        storedAbsFrames = expectedAbsFrames;
    end
    if isfield(SS, 'reconMeta') && ~isempty(SS.reconMeta)
        D.reconMeta = SS.reconMeta;
        if ~(SS.reconMeta.reconMode == "mvperframe" && SS.reconMeta.source == "raw")
            warning('mcLoadDataset:provenance', ...
                'storedOptimRaw.mat in %s reports reconMode=%s, source=%s (expected mvperframe/raw).', ...
                parent_folder, SS.reconMeta.reconMode, SS.reconMeta.source);
        end
    else
        D.reconMeta = struct('reconMode',"unknown",'source',"unknown");
    end

    absFrames = storedAbsFrames;
    startI    = absFrames(1);
    endI      = absFrames(end);
    nframes   = numel(absFrames);

    D.startI  = startI;
    D.endI    = endI;
    D.absFrames = absFrames;
    D.nframes = nframes;

    % -------------------- Load + undistort labels, build tensors --------------------
    labels0   = nan(nframes, nparts, ncams, 2);
    llh       = zeros(nframes, nparts, ncams);
    wantedLLH = {'l_eye_LLH','r_eye_LLH','nose_LLH','l_ear_LLH','r_ear_LLH', ...
                 'spine_LLH','neck_base_LLH'};

    for ii = 1:ncams
        cam       = camIDs(ii);
        camLabels = readtable(fullfile(parent_folder, sprintf('cam%d_labels_wide.csv', cam)));
        LLHtbl    = readtable(fullfile(parent_folder, sprintf('cam%d_llh_wide.csv',    cam)));

        if ismember('frame', camLabels.Properties.VariableNames)
            camLabels = removevars(camLabels, 'frame');
        end
        if ismember('frame', LLHtbl.Properties.VariableNames)
            LLHtbl = removevars(LLHtbl, 'frame');
        end

        row    = find(O.cam == cam, 1);
        assert(~isempty(row), 'Camera %d not found in offset.xlsx.', cam);
        offset = O.offset(row);

        rowStart = startI + offset;
        rowEnd   = endI   + offset;
        assert(rowStart >= 1 && rowEnd <= height(camLabels) && rowEnd <= height(LLHtbl), ...
            'cam%d crop range [%d %d] out of bounds.', cam, rowStart, rowEnd);

        camLabelsCrop = camLabels(rowStart:rowEnd, :);
        LLHcrop       = LLHtbl(rowStart:rowEnd, :);

        [camUDx, camUDy] = undistortBlenderLabels(camLabelsCrop, intrin{ii}, nframes, varnames);

        for p = 1:nparts
            labels0(:, p, ii, 1) = camUDx.(partNames{p});
            labels0(:, p, ii, 2) = camUDy.(partNames{p});
        end

        LLHsub = LLHcrop(:, wantedLLH);
        for p = 1:nparts
            llh(:, p, ii) = LLHsub.(wantedLLH{p});
        end
    end

    D.labels0   = labels0;
    D.llh       = llh;
    D.camIDs    = camIDs;
    D.partNames = partNames;
    D.nparts    = nparts;
    D.ncams     = ncams;

    % -------------------- Ear geometry (nominal offsets + measured scatter) --------------------
    % Nominal ear local offsets MUST match cheetah3D so the nominal recon
    % equals production. The per-axis scatter (used to perturb the ears in the
    % MC) is taken from earGeometry_results.mat (triLocal: the triangulated ear
    % position in the head-local frame across the few calibration frames).
    D.earOffL0 = [-0.0725,  0.0608, 0.0441];   % == plear in cheetah3D.m (production)
    D.earOffR0 = [-0.0725, -0.0608, 0.0441];   % == prear in cheetah3D.m (production)
    D.earStdL  = [];  D.earStdR = [];  D.earStdFrames = 0;
    egf = fullfile(parent_folder, 'earGeometry', 'earGeometry_results.mat');
    if isfile(egf)
        EG = load(egf);
        if isfield(EG, 'results') && isfield(EG.results, 'triLocal')
            tl = EG.results.triLocal;
            if isfield(tl, 'l_ear')
                Ll = tl.l_ear; Ll = Ll(all(isfinite(Ll),2), :);
                if size(Ll,1) >= 2, D.earStdL = std(Ll, 0, 1); D.earStdFrames = size(Ll,1); end
            end
            if isfield(tl, 'r_ear')
                Lr = tl.r_ear; Lr = Lr(all(isfinite(Lr),2), :);
                if size(Lr,1) >= 2, D.earStdR = std(Lr, 0, 1); end
            end
        end
    end

    % -------------------- Lure 3D positions (for alpha/beta/gamma/delta) --------------------
    % Align lure samples to the absolute reconstruction frames.
    % (intersect cheetah reconFrames with Lure3D.Index). Held fixed at the
    % nominal triangulation; lure-position uncertainty is a separate, smaller
    % term (the lure is metres away, so it perturbs the lure DIRECTION little).
    lureXYZ = nan(nframes, 3);
    lf = dir(fullfile(parent_folder, '*LurePos.mat'));
    if ~isempty(lf)
        LL = load(fullfile(parent_folder, lf(1).name));
        if isfield(LL, 'Lure3D')
            l3 = LL.Lure3D;
            [~, rLoc, lLoc] = intersect(absFrames, l3.Index(:), 'stable');
            marks = l3.Mark(lLoc);
            for q = 1:numel(rLoc)
                v = marks{q};
                lureXYZ(rLoc(q), :) = v(:).';
            end
        else
            warning('%s has no Lure3D field; alpha/beta/gamma/delta will be NaN.', lf(1).name);
        end
    else
        warning('No *LurePos.mat found; alpha/beta/gamma/delta will be NaN.');
    end
    D.lureXYZ   = lureXYZ;
    D.lureValid = all(isfinite(lureXYZ), 2);

    % -------------------- Per-camera RAW lure 2D labels (for lure perturbation) --------------------
    % cam<ID>lureLabels.mat stores the RAW (distorted) lure pixels + absolute
    % Index, gated identically to production (lureReconTriangulation.m). We load
    % them so the MC can perturb the lure labels in raw px, undistort, and
    % re-triangulate (mcLureTriangulate) -- the lure-position uncertainty source.
    % lureLabels2D is NaN where the lure was not visible in that camera/frame.
    lureLabels2D = nan(nframes, 2, ncams);
    for ii = 1:ncams
        cam = camIDs(ii);
        llFile = fullfile(parent_folder, sprintf('cam%dlureLabels.mat', cam));
        if ~isfile(llFile)
            warning('mcLoadDataset:noLureLabels', ...
                'cam%dlureLabels.mat not found; that camera will not contribute to lure perturbation.', cam);
            continue
        end
        LR = load(llFile, 'results');
        if ~isfield(LR, 'results') || isempty(LR.results), continue; end
        res = LR.results;
        rows = res.Index(:) - startI + 1;            % absolute Index -> recon row
        keep = rows >= 1 & rows <= nframes;
        rows = rows(keep);
        lureLabels2D(rows, 1, ii) = res.X(keep);
        lureLabels2D(rows, 2, ii) = res.Y(keep);
    end
    D.lureLabels2D = lureLabels2D;                    % [nframes x 2 x ncams] raw distorted XY
    D.lureVis      = squeeze(all(isfinite(lureLabels2D), 2));  % [nframes x ncams] visibility
    if ncams == 1, D.lureVis = D.lureVis(:); end     % keep [nframes x 1] shape

    % -------------------- Per-part pixel sigma (matches reconPerFrame.m) --------------------
    %   l_eye r_eye nose l_ear r_ear spine neck_base
    D.sigmas = [3 3 2 5 5 6 6];

    % -------------------- Per-frame camera-trust mask (from multiview LOO pass) --------
    % Port of reconPerFrame.m:326-350. multiviewRecon.m flagged cameras carrying a
    % systematic bias per frame via leave-one-out reprojection (camKeptMask, rows
    % aligned to absFrames). A false entry HARD-EXCLUDES that camera from the frame's
    % fit (same as LLH=0), so the MC uses the identical camera set as production.
    % Missing file/field => keep all cameras (graceful fallback, warn once).
    camMaskStore = true(nframes, ncams);
    camLOOStore  = nan(nframes, ncams);
    triFile = fullfile(parent_folder, 'triangulatedPointsRaw_multiview.mat');
    if isfile(triFile)
        TRI = load(triFile, 'camKeptMask', 'camLOOResid', 'camIDs', 'absFrames');
        if isfield(TRI, 'camKeptMask') && isfield(TRI, 'camIDs') && isfield(TRI, 'absFrames')
            triAbs   = double(TRI.absFrames(:));
            triRowOf = containers.Map('KeyType','double','ValueType','double');
            for rr = 1:numel(triAbs)
                triRowOf(triAbs(rr)) = rr;
            end
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
            fprintf('Per-frame camera-trust mask loaded: %d camera-frames rejected (of %d).\n', ...
                sum(~camMaskStore(:)), numel(camMaskStore));
        else
            warning('mcLoadDataset:noCamMask', ...
                '%s has no camKeptMask/camIDs/absFrames; using all cameras in the fit.', triFile);
        end
    else
        warning('mcLoadDataset:noTriFile', ...
            'triangulatedPointsRaw_multiview.mat not found in %s; using all cameras in the fit.', parent_folder);
    end
    D.camMaskStore = camMaskStore;
    D.camLOOStore  = camLOOStore;
end
