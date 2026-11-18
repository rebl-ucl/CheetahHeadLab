function R = mcMainAnalysis(datasetName, opts)
% MCMAINANALYSIS  Clean main Monte-Carlo angle-uncertainty analysis (CRBR).
%
%   R = mcMainAnalysis(datasetName, opts)
%
% Reusable across datasets: pass the dataset folder name (under
% <repo>/outputs by default) and it runs the CRBR reconstruction
% Monte-Carlo, reporting the per-frame uncertainty of the head-stabilisation
% angles. The primary/default estimator is reconTrajectoryBatch: every draw
% receives a robust per-frame initial fit followed by the same bidirectional
% SO(3) trajectory reconstruction used by the production pipeline.
%
% THREE uncertainty sources (kept deliberately explicit):
%   1. 2D LABEL NOISE   - iid Gaussian, configured by
%      opts.perturbRadialPx, per frame/point.
%   2. KINEMATIC LENGTHS - every head+torso template length is a variable with
%      a NOMINAL value (= the reconstruction value, template x estimated scale)
%      and a per-length 1-sigma. One draw per MC sample is shared across all
%      frames (the animal's geometry is fixed but uncertain). Edit the sigmas
%      via opts.lengthSigma / opts.lengthSigmaFrac.
%   3. LURE LABEL NOISE - visible 2-D lure labels are perturbed and the 3-D
%      lure is re-triangulated independently in every draw.
%
% Each MC sample first re-fits every frame's 12-DOF pose independently with
% the production reprojection cost + robust loss, then reconstructs the full
% trajectory bidirectionally. It records six angles (mcMainAngles):
% gamma/delta (head/body yaw error),
% alpha/beta (head/body pitch error), and headLureAngle/bodyLureAngle (the
% TOTAL gaze-to-lure angle, a single off-target sanity check). Aggregated to
% per-frame median + 95% CI.
%
% opts (all optional):
%   dataRoot        <repo>/outputs
%   nMC             500
%   perturbRadialPx [2.5 2.5 2.5 3 3 3.5 3.5]  per-part label-noise RMS on the TOTAL
%                            2D displacement (px), order partNames
%                            {l_eye,r_eye,nose,l_ear,r_ear,spine,neck_base}.
%                            Per-axis sigma = perturbRadialPx/sqrt(2).
%   perturbScale    1.0      global multiplier on perturbRadialPx (sweeps)
%   lengthSigmaFrac 0.08     default per-length 1-sigma as a fraction of |nominal|
%                            (fallback for any length not in lengthSigma)
%   lengthSigma     struct() per-length sigma override in METRES, e.g.
%                            struct('ear_x',0.005,'spine_x',0.01)
%   clampLengthSign false    keep each perturbed length the sign of its nominal
%   lurePerturbPx   2        per-axis lure-label sigma; 0 conditions on nominal lure
%   camIDs          1:6
%   robustType/robustC  'huber'/3  (consensus rejection auto-ON when robust~='none',
%                            matching production cheetahError: C=4, floorSigma=2, minCams=3)
%   reconstructionMode 'trajectory' (primary); 'perframe' is retained for audit
%   trajectoryNOuter 6      bidirectional reconstruction passes per draw
%   despikeAngles   false    optional Hampel + Savitzky-Golay post-processing;
%                            raw angles are always retained independently
%   hampelK/hampelSigma  3/3  Hampel half-window / MAD threshold
%   sgolayOrder/sgolayWin 3/9  Savitzky-Golay order / frame length
%   useParallel     true
%   seed            1
%   makePlots       true
%   outDir          <dataset>/MonteCarloResults
%   tag             'main'

    % ---------- defaults ----------
    if nargin < 2, opts = struct(); end
    thisDir  = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(thisDir));
    def = struct('dataRoot', fullfile(repoRoot,'outputs'), ...
        'nMC',500,'lengthSigmaFrac',0.08, ...
        'perturbRadialPx',[2.5 2.5 2.5 3 3 3.5 3.5],'perturbScale',1, ...
        'clampLengthSign',false, 'lurePerturbPx',2, ...
        'camIDs',1:6,'robustType','huber','robustC',3, ...
        'reconstructionMode','trajectory','trajectoryNOuter',6, ...
        'trajectoryAngleMidSigmaDeg',2.0,'trajectoryPositionMidSigmaM',0.015, ...
        'trajectoryTemporalHuberC',3.0, ...
        'despikeAngles',false,'hampelK',3,'hampelSigma',3,'sgolayOrder',3,'sgolayWin',9, ...
        'savePoseSamples',true,'savePerFramePoseSamples',true, ...
        'saveLureSamples',true,'saveFitDiagnostics',true, ...
        'useParallel',true,'seed',1, ...
        'makePlots',true,'outDir','','tag','main');
    fnD = fieldnames(def);
    for k = 1:numel(fnD)
        if ~isfield(opts, fnD{k}) || isempty(opts.(fnD{k})), opts.(fnD{k}) = def.(fnD{k}); end
    end
    if ~isfield(opts,'lengthSigma')
        opts.lengthSigma = struct('ear_x',0.011,'ear_y',0.005,'ear_z',0.009);
    elseif ~isstruct(opts.lengthSigma)
        error('lengthSigma must be a struct of absolute length sigmas in metres.');
    end
    addpath(fullfile(repoRoot,'CRBR','reconstruction'));
    parent_folder = fullfile(opts.dataRoot, datasetName);
    assert(isfolder(parent_folder), 'Dataset folder not found: %s', parent_folder);
    if isempty(opts.outDir), opts.outDir = fullfile(parent_folder,'MonteCarloResults'); end
    if ~exist(opts.outDir,'dir'), mkdir(opts.outDir); end

    Rx = @(t)[1 0 0; 0 cos(t) sin(t); 0 -sin(t) cos(t)];
    Ry = @(t)[cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
    Rz = @(t)[cos(t) sin(t) 0; -sin(t) cos(t) 0; 0 0 1];
    % Robust loss + per-landmark consensus rejection, matching production
    % (cheetahError.m): consensus is ON whenever a robust loss is active.
    robustOpt = struct('type', opts.robustType, 'c', opts.robustC, ...
        'consensus', struct('enable', ~strcmp(opts.robustType,'none'), ...
                            'C',4, 'floorSigma',2, 'minCams',3));

    angleNames = {'alpha_headPitchErr','beta_bodyPitchErr', ...
                  'gamma_headYawErr','delta_bodyYawErr', ...
                  'headLureAngle','bodyLureAngle'};
    nAng = numel(angleNames);

    % ---------- load dataset ----------
    fprintf('Dataset: %s\n', parent_folder);
    D = mcLoadDataset(parent_folder, opts.camIDs);
    nframes = D.nframes; ncams = D.ncams; nparts = D.nparts;
    partNames = D.partNames; sigmas = D.sigmas;
    labels0 = D.labels0; llh = D.llh; lureXYZ = D.lureXYZ;
    invPoses = cellfun(@(P) invert(P), D.poses, 'UniformOutput', false);
    intrin = D.intrin;
    camMaskStore = D.camMaskStore;              % [nframes x ncams] LOO camera-trust mask
    absFrames = D.absFrames(:);

    % ---------- lure-perturbation inputs (third uncertainty source) ----------
    % Built once; identical to lureReconTriangulation.m's cameraPoses table and
    % intrinsics array. When opts.lurePerturbPx==0 these are unused (lure held
    % fixed at D.lureXYZ, i.e. the original behaviour).
    lurePerturbPx = opts.lurePerturbPx;
    lureLabels2D  = D.lureLabels2D;            % [nframes x 2 x ncams] raw distorted XY
    lureVis       = D.lureVis;                 % [nframes x ncams] visibility
    camPosesTbl   = table(uint32(D.camIDs(:)), vertcat(D.poses{:}), ...
                          'VariableNames', {'ViewId','AbsolutePose'});
    intrinArr     = [D.intrin{:}];
    lureCamIDs    = D.camIDs;

    % ---------- kinematic length table ----------
    % Template offsets MUST match cheetah3D (the production reconstruction cost
    % cheetahError->cheetah3D), so the nominal MC fit reproduces production. The
    % ear offsets below equal cheetah3D's plear/prear [-0.0725 0.0608 0.0441];
    % they were previously [-0.0707 0.0457 0.0479], a ~1.5 cm (33%) mismatch.
    lengthNames = {'eye_y','nose_x','nose_z','ear_x','ear_y','ear_z', ...
                   'neck_x','neck_z','spine_x','spine_z'};
    template = [0.03, 0.055, -0.055, -0.0725, 0.0608, 0.0441, ...
                0.04, 0.10, -0.33, 0.10];
    isHead   = logical([1 1 1 1 1 1 0 0 0 0]);
    scaleVec = isHead*D.scaleHead0 + ~isHead*D.scaleTorso0;
    Lnom     = template .* scaleVec;                      % nominal lengths (m)
    Lsig     = opts.lengthSigmaFrac * abs(Lnom);          % default sigmas (m)
    for k = 1:numel(lengthNames)                          % per-name overrides
        if isfield(opts.lengthSigma, lengthNames{k})
            Lsig(k) = opts.lengthSigma.(lengthNames{k});
        end
    end
    offNom = offsetsFromLengths(Lnom);

    fprintf('\nKinematic length table (nominal = template x estimated scale):\n');
    fprintf('  scaleHead=%.4f  scaleTorso=%.4f\n', D.scaleHead0, D.scaleTorso0);
    fprintf('  %-9s %12s %12s %10s\n','length','nominal(m)','sigma(m)','sigma(%)');
    for k = 1:numel(lengthNames)
        fprintf('  %-9s %12.4f %12.4f %9.1f%%\n', lengthNames{k}, Lnom(k), Lsig(k), ...
            100*Lsig(k)/max(abs(Lnom(k)),eps));
    end
    % Per-part label-noise sigma. perturbRadialPx is the sigma on the TOTAL 2D
    % displacement; per-axis sigma = radial/sqrt(2) (independent x and y). Fit
    % weights stay = D.sigmas (production variances), unchanged.
    perAxisSig = (opts.perturbScale * opts.perturbRadialPx(:)) / sqrt(2);
    assert(numel(perAxisSig)==nparts, 'perturbRadialPx must have %d entries (partNames order).', nparts);
    fprintf('\nLabel perturbation (radial px, per part): %s\n', mat2str(opts.perturbScale*opts.perturbRadialPx,3));
    fprintf('Frames: %d (abs %d..%d) | nMC=%d | lure for %d frames\n', ...
        nframes, absFrames(1), absFrames(end), opts.nMC, nnz(D.lureValid));
    if lurePerturbPx > 0
        fprintf('Lure perturbation: %.1f px (isotropic, per-frame iid) + re-triangulation\n', lurePerturbPx);
    else
        fprintf('Lure perturbation: OFF (held fixed at nominal triangulation)\n');
    end

    % ---------- solver setup ----------
    % Match the production per-frame main solve (reconPerFrame.m:530-535): loose
    % tolerances with a high iteration cap (the cap is a safety net; lsqnonlin
    % converges well before it), so the MC solves to the same criterion.
    solveOpts = optimoptions(@lsqnonlin,'Display','off', ...
        'StepTolerance',1e-4,'OptimalityTolerance',1e-4, ...
        'MaxFunctionEvaluations',1e5,'MaxIterations',1e5);

    % Production physiological bounds (reconPerFrame.m:454-472): asymmetric angle
    % limits + camera-extent position box. GLOBAL bounds (not a local box around
    % the seed), identical to the reconstruction the MC characterises.
    xcams = cellfun(@(P) P.Translation(1), D.poses);
    ycams = cellfun(@(P) P.Translation(2), D.poses);
    zcams = cellfun(@(P) P.Translation(3), D.poses);
    LB = zeros(1,12); UB = zeros(1,12);
    LB(1)  = -pi/5;            UB(1)  =  pi/5;
    LB(2)  = -pi/3;            UB(2)  =  pi/3;
    LB(3)  = -pi;              UB(3)  =  pi;
    LB(4)  = min(xcams) - 5;   UB(4)  = max(xcams) + 5;
    LB(5)  = min(ycams) - 2;   UB(5)  = max(ycams) + 2;
    LB(6)  = min(zcams) - 0.5; UB(6)  = max(zcams) + 1.5;
    LB(7)  = -pi/5;            UB(7)  =  pi/5;
    LB(8)  = -pi/3;            UB(8)  =  pi/3;
    LB(9)  = -pi;              UB(9)  =  pi;
    LB(10) = min(xcams) - 5;   UB(10) = max(xcams) + 5;
    LB(11) = min(ycams) - 2;   UB(11) = max(ycams) + 2;
    LB(12) = min(zcams) - 0.5; UB(12) = max(zcams) + 1.5;

    % Optional despiking settings for the angle series:
    % the per-frame fit has no temporal consistency, so weak-evidence frames make
    % isolated spikes in any derived angle. Despike all six outputs per MC
    % sample and for the nominal so component and total-angle panels have
    % consistent temporal processing. Total angles are unsigned and remain in
    % [0,pi] after smoothing; the four signed errors remain wrapped to +/-pi.
    dsp = struct('on', opts.despikeAngles, 'K', opts.hampelK, 'nSigma', opts.hampelSigma, ...
                 'sgOrder', opts.sgolayOrder, 'sgWin', opts.sgolayWin, ...
                 'cols', 1:nAng, 'unsignedCols', 5:6);

    reconstructionMode = lower(string(opts.reconstructionMode));
    if ismember(reconstructionMode, ["trajectory","bidirectional","robust_bidirectional"])
        reconstructionMode = "trajectory";
    elseif reconstructionMode ~= "perframe"
        error('Unknown reconstructionMode "%s". Use trajectory or perframe.', reconstructionMode);
    end

    % ---------- nominal pass ----------
    % The nominal IS the production reconstruction: the reconPerFrame poses in
    % storedOptimRaw.mat and, for the primary estimator, the reconTrajectoryBatch
    % poses in storedOptimRaw_trajectory.mat. Re-estimating it would not
    % reproduce the pipeline exactly: the per-frame problem is flat near its
    % minimum, so re-solving production's own cost from its own solution moves
    % a median ~0.3 deg. The nominal trajectory is still re-run through this code
    % path (the one every draw uses) and its difference is reported.
    fprintf('\nNominal = production reconstruction (storedOptimRaw.mat');
    assert(size(D.xBarStore,1) == nframes, 'storedOptimRaw.mat rows do not match the frame range.');
    nominalPerFramePose = D.xBarStore;
    nominalPerFrameRMSE = poseFitRMSE(nominalPerFramePose, D, offNom, Rz,Rx,Ry, ...
        sigmas, invPoses, intrin, ncams, nparts, partNames, robustOpt);

    nominalTrajectoryDiagnostics = table();
    nominalRerunPose = [];
    nominalRerunMaxDiffDeg = NaN;
    if reconstructionMode == "trajectory"
        fprintf(', storedOptimRaw_trajectory.mat)\n');
        trajFile = fullfile(parent_folder, 'storedOptimRaw_trajectory.mat');
        assert(isfile(trajFile), '%s is missing: run reconTrajectoryBatch (runFullPipeline) first.', trajFile);
        TP = load(trajFile, 'xBarStore', 'absFrames');
        assert(isequal(double(TP.absFrames(:)), absFrames), ...
            '%s frames do not match storedOptimRaw.mat: re-run the pipeline.', trajFile);
        nominalPose = TP.xBarStore;
        fprintf('Reproducibility check: nominal trajectory re-run through the draw code path...\n');
        tOpts = makeTrajectoryOpts(opts, D, nominalPerFramePose, offNom);
        Tnom = reconTrajectoryBatch(parent_folder, tOpts);
        nominalRerunPose = Tnom.pose;
        nominalTrajectoryDiagnostics = Tnom.diagnostics;
        angCols = [1 2 3 7 8 9];
        nominalRerunMaxDiffDeg = rad2deg(max(abs(wrapToPi( ...
            nominalRerunPose(:,angCols) - nominalPose(:,angCols))), [], 'all'));
        fprintf('  re-run vs production trajectory: max %.4f deg\n', nominalRerunMaxDiffDeg);
    else
        fprintf(')\n');
        nominalPose = nominalPerFramePose;
    end

    nominalAnglesRaw = anglesForTrajectory(nominalPose, lureXYZ, Rz,Rx,Ry, nAng);
    nominalAnglesProcessed = processAngles(nominalAnglesRaw, dsp);
    nominalTrajectoryRMSE = poseFitRMSE(nominalPose, D, offNom, Rz,Rx,Ry, ...
        sigmas, invPoses, intrin, ncams, nparts, partNames, robustOpt);
    nominalPerFrameAnglesRaw = anglesForTrajectory(nominalPerFramePose, lureXYZ, Rz,Rx,Ry, nAng);
    fprintf(['Nominal per-frame/primary reproj RMSE (sigma units): median ' ...
        '%.3f / %.3f; p95 %.3f / %.3f\n'], ...
        median(nominalPerFrameRMSE,'omitnan'), median(nominalTrajectoryRMSE,'omitnan'), ...
        prctile(nominalPerFrameRMSE,95), prctile(nominalTrajectoryRMSE,95));

    % ---------- Monte-Carlo ----------
    nMC = opts.nMC; clampSign = opts.clampLengthSign;
    rng(opts.seed,'twister');
    streams = RandStream.create('mrg32k3a','NumStreams',nMC,'CellOutput',true,'Seed',opts.seed);
    anglesRawAll = nan(nMC, nframes, nAng, 'single');
    anglesProcessedAll = nan(nMC, nframes, nAng, 'single');
    lengthDraws = nan(nMC, numel(Lnom));
    sampleSucceeded = false(nMC,1);
    sampleError = cell(nMC,1);
    if opts.savePoseSamples
        poseAll = nan(nMC, nframes, 12, 'single');
    else
        poseAll = single([]);
    end
    if opts.savePerFramePoseSamples
        perFramePoseAll = nan(nMC, nframes, 12, 'single');
    else
        perFramePoseAll = single([]);
    end
    if opts.saveLureSamples
        lureXYZAll = nan(nMC, nframes, 3, 'single');
    else
        lureXYZAll = single([]);
    end
    if opts.saveFitDiagnostics
        perFrameFitRMSEAll = nan(nMC,nframes,'single');
        trajectoryFitRMSEAll = nan(nMC,nframes,'single');
    else
        perFrameFitRMSEAll = single([]);
        trajectoryFitRMSEAll = single([]);
    end
    trajectoryIterations = zeros(nMC,1,'uint16');
    trajectoryFinalMedianRotationUpdateDeg = nan(nMC,1);
    trajectoryFinalP95RotationUpdateDeg = nan(nMC,1);
    trajectoryFinalMedianPositionUpdateM = nan(nMC,1);

    usedPar = false;
    if opts.useParallel && ~isempty(ver('parallel'))
        try
            if isempty(gcp('nocreate')), parpool; end
            usedPar = true;
        catch ME
            warning('No parallel pool (%s); running serially.', ME.message);
        end
    end

    lureArgs = {lurePerturbPx, lureLabels2D, lureVis, camPosesTbl, intrinArr, lureCamIDs};
    tS = tic;
    if usedPar
        parfor mc = 1:nMC
            S = oneSample(streams{mc}, D, parent_folder, nominalPose, labels0, llh, ...
                perAxisSig, clampSign, Lnom, Lsig, Rz,Rx,Ry, sigmas, invPoses, intrin, ...
                ncams, nparts, partNames, robustOpt, LB, UB, camMaskStore, solveOpts, ...
                nAng, lureXYZ, lureArgs, dsp, opts, reconstructionMode);
            anglesRawAll(mc,:,:) = reshape(single(S.anglesRaw),[1 nframes nAng]);
            anglesProcessedAll(mc,:,:) = reshape(single(S.anglesProcessed),[1 nframes nAng]);
            lengthDraws(mc,:) = S.lengthDraw;
            sampleSucceeded(mc) = S.succeeded;
            sampleError{mc} = S.error;
            if opts.savePoseSamples, poseAll(mc,:,:) = reshape(single(S.pose),[1 nframes 12]); end
            if opts.savePerFramePoseSamples, perFramePoseAll(mc,:,:) = reshape(single(S.perFramePose),[1 nframes 12]); end
            if opts.saveLureSamples, lureXYZAll(mc,:,:) = reshape(single(S.lureXYZ),[1 nframes 3]); end
            if opts.saveFitDiagnostics
                perFrameFitRMSEAll(mc,:) = single(S.perFrameFitRMSE);
                trajectoryFitRMSEAll(mc,:) = single(S.trajectoryFitRMSE);
            end
            trajectoryIterations(mc) = uint16(S.trajectoryIterations);
            trajectoryFinalMedianRotationUpdateDeg(mc) = S.finalMedianRotationUpdateDeg;
            trajectoryFinalP95RotationUpdateDeg(mc) = S.finalP95RotationUpdateDeg;
            trajectoryFinalMedianPositionUpdateM(mc) = S.finalMedianPositionUpdateM;
        end
    else
        for mc = 1:nMC
            S = oneSample(streams{mc}, D, parent_folder, nominalPose, labels0, llh, ...
                perAxisSig, clampSign, Lnom, Lsig, Rz,Rx,Ry, sigmas, invPoses, intrin, ...
                ncams, nparts, partNames, robustOpt, LB, UB, camMaskStore, solveOpts, ...
                nAng, lureXYZ, lureArgs, dsp, opts, reconstructionMode);
            anglesRawAll(mc,:,:) = reshape(single(S.anglesRaw),[1 nframes nAng]);
            anglesProcessedAll(mc,:,:) = reshape(single(S.anglesProcessed),[1 nframes nAng]);
            lengthDraws(mc,:) = S.lengthDraw;
            sampleSucceeded(mc) = S.succeeded;
            sampleError{mc} = S.error;
            if opts.savePoseSamples, poseAll(mc,:,:) = reshape(single(S.pose),[1 nframes 12]); end
            if opts.savePerFramePoseSamples, perFramePoseAll(mc,:,:) = reshape(single(S.perFramePose),[1 nframes 12]); end
            if opts.saveLureSamples, lureXYZAll(mc,:,:) = reshape(single(S.lureXYZ),[1 nframes 3]); end
            if opts.saveFitDiagnostics
                perFrameFitRMSEAll(mc,:) = single(S.perFrameFitRMSE);
                trajectoryFitRMSEAll(mc,:) = single(S.trajectoryFitRMSE);
            end
            trajectoryIterations(mc) = uint16(S.trajectoryIterations);
            trajectoryFinalMedianRotationUpdateDeg(mc) = S.finalMedianRotationUpdateDeg;
            trajectoryFinalP95RotationUpdateDeg(mc) = S.finalP95RotationUpdateDeg;
            trajectoryFinalMedianPositionUpdateM(mc) = S.finalMedianPositionUpdateM;
            if mod(mc,max(1,round(nMC/10)))==0, fprintf('  MC %d/%d (%.0fs)\n', mc,nMC,toc(tS)); end
        end
    end
    fprintf('Monte-Carlo done in %.1f s (parallel=%d; successful=%d/%d).\n', ...
        toc(tS), usedPar, nnz(sampleSucceeded), nMC);

    % ---------- aggregate (wrap-safe, deviations from nominal) ----------
    processedSummary = aggregateAngles(double(anglesProcessedAll), nominalAnglesProcessed);
    rawSummary = aggregateAngles(double(anglesRawAll), nominalAnglesRaw);

    R = struct();
    R.datasetName = datasetName; R.parent_folder = parent_folder; R.opts = opts;
    % Nominal reconstruction provenance: which trajectory the MC perturbs around
    % (from storedOptimRaw.mat's reconMeta) and the robust loss used here.
    estimatorName = "robust_bidirectional_trajectory";
    if reconstructionMode == "perframe", estimatorName = "robust_independent_perframe"; end
    R.provenance = struct('reconMeta', D.reconMeta, ...
        'estimator', estimatorName, ...
        'reconstructionMode', reconstructionMode, ...
        'robustType', string(opts.robustType), 'robustC', opts.robustC, ...
        'consensus', robustOpt.consensus, ...
        'camMaskSource', "triangulatedPointsRaw_multiview.mat:camKeptMask", ...
        'camFramesRejected', sum(~camMaskStore(:)), ...
        'despike', struct('on', dsp.on, 'hampelK', dsp.K, 'hampelSigma', dsp.nSigma, ...
                          'sgolayOrder', dsp.sgOrder, 'sgolayWin', dsp.sgWin, ...
                          'cols', dsp.cols, 'unsignedCols', dsp.unsignedCols));
    R.provenance.randomisation = struct('masterSeed',opts.seed, ...
        'generator',"mrg32k3a", 'sampleIndexDefinesStream',true, ...
        'labelNoiseModel',"independent Gaussian per axis", ...
        'lengthModel',"one Gaussian geometry draw per animal trajectory", ...
        'lureNoiseModel',"independent Gaussian per visible 2-D label axis");
    R.angleNames = angleNames; R.absFrames = absFrames; R.nframes = nframes;
    R.lengthNames = lengthNames; R.Lnom = Lnom; R.Lsig = Lsig;
    R.sampleIndex = (1:nMC).';
    R.camIDs = D.camIDs;
    R.camMaskStore = camMaskStore;
    R.partNames = partNames;
    R.labelPerAxisSigmaPx = perAxisSig;
    R.nominalOffsets = offNom;
    R.lurePerturbPx = lurePerturbPx;
    R.nominalDeg = processedSummary.nominalDeg;
    R.medianDeg = processedSummary.medianDeg;
    R.ciLoDeg = processedSummary.ciLoDeg;
    R.ciHiDeg = processedSummary.ciHiDeg;
    R.ciHalfDeg = processedSummary.ciHalfDeg;
    R.stdDeg = processedSummary.stdDeg;
    R.rawSummary = rawSummary;
    R.nominalAnglesRawDeg = rad2deg(nominalAnglesRaw);
    R.nominalAnglesProcessedDeg = rad2deg(nominalAnglesProcessed);
    R.nominalPerFrameAnglesRawDeg = rad2deg(nominalPerFrameAnglesRaw);
    R.nominalPose = nominalPose;
    R.nominalPerFramePose = nominalPerFramePose;
    R.nominalSource = "production storedOptimRaw.mat / storedOptimRaw_trajectory.mat";
    R.nominalRerunPose = nominalRerunPose;
    R.nominalRerunMaxDiffDeg = nominalRerunMaxDiffDeg;
    R.nominalRMSE = nominalTrajectoryRMSE;
    R.nominalPerFrameRMSE = nominalPerFrameRMSE;
    R.nominalTrajectoryRMSE = nominalTrajectoryRMSE;
    R.nominalLureXYZ = lureXYZ;
    R.nominalTrajectoryDiagnostics = nominalTrajectoryDiagnostics;
    R.anglesRawAllDeg = rad2deg(anglesRawAll);
    R.anglesProcessedAllDeg = rad2deg(anglesProcessedAll);
    R.anglesAllDeg = R.anglesProcessedAllDeg; % backwards-compatible alias
    R.poseAll = poseAll;
    R.perFramePoseAll = perFramePoseAll;
    R.lureXYZAll = lureXYZAll;
    R.lengthDraws = lengthDraws;
    R.perFrameFitRMSEAll = perFrameFitRMSEAll;
    R.trajectoryFitRMSEAll = trajectoryFitRMSEAll;
    R.sampleSucceeded = sampleSucceeded;
    R.sampleError = sampleError;
    R.trajectoryIterations = trajectoryIterations;
    R.trajectoryFinalMedianRotationUpdateDeg = trajectoryFinalMedianRotationUpdateDeg;
    R.trajectoryFinalP95RotationUpdateDeg = trajectoryFinalP95RotationUpdateDeg;
    R.trajectoryFinalMedianPositionUpdateM = trajectoryFinalMedianPositionUpdateM;

    % ---------- report ----------
    fprintf('\n=== Per-frame angle uncertainty (95%% CI half-width, deg) ===\n');
    fprintf('%-20s  median   p90    max\n','angle');
    for a = 1:nAng
        c = R.ciHalfDeg(:,a);
        fprintf('%-20s  %5.2f   %5.2f  %5.2f\n', angleNames{a}, ...
            median(c,'omitnan'), prctile(c,90), max(c));
    end

    % ---------- save ----------
    tag = opts.tag; if ~isempty(tag) && tag(1)~='_', tag=['_' tag]; end
    save(fullfile(opts.outDir, sprintf('mcMain%s.mat', tag)), '-struct','R','-v7.3');
    writeCsv(R, fullfile(opts.outDir, sprintf('mcMain_per_frame%s.csv', tag)));
    if opts.makePlots, makeMainFigures(R, opts.outDir, tag); end
    fprintf('\nOutputs in: %s\n', opts.outDir);
end

% ========================================================================
function S = oneSample(stream, D, parentFolder, nominalPose, labels0, llh, perAxisSig, ...
        clampSign, Lnom, Lsig, Rz,Rx,Ry, sigmas, invPoses, intrin, ncams, nparts, ...
        partNames, robustOpt, LB, UB, camMaskStore, solveOpts, nAng, lureXYZ, lureArgs, ...
        dsp, opts, reconstructionMode)
    nframes = D.nframes;
    S = struct('anglesRaw',nan(nframes,nAng), 'anglesProcessed',nan(nframes,nAng), ...
        'pose',nan(nframes,12), 'perFramePose',nan(nframes,12), ...
        'lureXYZ',nan(nframes,3), 'lengthDraw',nan(1,numel(Lnom)), ...
        'perFrameFitRMSE',nan(1,nframes), 'trajectoryFitRMSE',nan(1,nframes), ...
        'trajectoryIterations',0, 'finalMedianRotationUpdateDeg',nan, ...
        'finalP95RotationUpdateDeg',nan, 'finalMedianPositionUpdateM',nan, ...
        'succeeded',false, 'error','');

    % One geometry draw is shared across every frame of this animal trajectory.
    Ldraw = Lnom + Lsig .* randn(stream,1,numel(Lnom));
    if clampSign
        pos = Lnom >= 0;
        Ldraw(pos) = max(Ldraw(pos),0);
        Ldraw(~pos) = min(Ldraw(~pos),0);
    end
    S.lengthDraw = Ldraw;
    off = offsetsFromLengths(Ldraw);

    [lurePerturbPx, lureLabels2D, lureVis, camPosesTbl, intrinArr, lureCamIDs] = lureArgs{:};
    try
        % Draw the whole observation tensor first: the temporal reconstruction
        % therefore sees one coherent noisy sequence, not independently filtered
        % frame outputs.
        labelNoise = randn(stream,[nframes nparts ncams 2]) .* ...
            reshape(perAxisSig,[1 nparts 1 1]);
        labelsDraw = labels0 + labelNoise;

        % Per-frame stage = production reconPerFrame structure: fresh seed (the
        % production per-frame solution) + this draw's own previous-frame
        % solution, keep the lower robust cost.
        perFramePose = nan(nframes,12);
        perFrameRMSE = nan(1,nframes);
        nRows = size(D.xBarStore,1);
        for f = 1:nframes
            seed = min(max(D.xBarStore(min(f,nRows),:), LB), UB);
            labF = squeeze(labelsDraw(f,:,:,:));
            llhF = squeeze(llh(f,:,:));
            cmask = camMaskStore(f,:);
            fobj = @(pp) mcMainResidual(pp, Rz,Rx,Ry, off, labF, llhF, ...
                sigmas, invPoses, intrin, ncams, nparts, partNames, robustOpt, cmask);
            if f == 1
                starts = {seed};
            else
                starts = {seed, min(max(perFramePose(f-1,:), LB), UB)};
            end
            psol = seed; bestRes = inf;
            for sIdx = 1:numel(starts)
                [cand,resn] = lsqnonlin(fobj,starts{sIdx},LB,UB,solveOpts);
                if resn < bestRes, bestRes = resn; psol = cand; end
            end
            perFramePose(f,:) = psol;
            perFrameRMSE(f) = sqrt(bestRes / max(2*nnz(llhF>0.5),1));
        end

        lureDraw = lureXYZ;
        if lurePerturbPx > 0
            for f = 1:nframes
                lab2D = reshape(lureLabels2D(f,:,:),2,[])';
                visF = logical(lureVis(f,:));
                lab2Dp = lab2D + lurePerturbPx*randn(stream,size(lab2D));
                lureDraw(f,:) = mcLureTriangulate(lab2Dp,visF(:), ...
                    camPosesTbl,intrinArr,lureCamIDs);
            end
        end

        Ddraw = D;
        Ddraw.labels0 = labelsDraw;
        if reconstructionMode == "trajectory"
            tOpts = makeTrajectoryOpts(opts,Ddraw,perFramePose,off);
            T = reconTrajectoryBatch(parentFolder,tOpts);
            pose = T.pose;
            diag = T.diagnostics;
            S.trajectoryIterations = height(diag);
            if ~isempty(diag)
                S.finalMedianRotationUpdateDeg = diag.median_rotation_update_deg(end);
                S.finalP95RotationUpdateDeg = diag.p95_rotation_update_deg(end);
                S.finalMedianPositionUpdateM = diag.median_position_update_m(end);
            end
        else
            pose = perFramePose;
        end

        S.perFramePose = perFramePose;
        S.pose = pose;
        S.lureXYZ = lureDraw;
        S.perFrameFitRMSE = perFrameRMSE;
        S.trajectoryFitRMSE = poseFitRMSE(pose,Ddraw,off,Rz,Rx,Ry, ...
            sigmas,invPoses,intrin,ncams,nparts,partNames,robustOpt).';
        S.anglesRaw = anglesForTrajectory(pose,lureDraw,Rz,Rx,Ry,nAng);
        S.anglesProcessed = processAngles(S.anglesRaw,dsp);
        S.succeeded = all(isfinite(pose(:))) && any(isfinite(S.anglesRaw(:)));
        if ~S.succeeded, S.error = 'Non-finite trajectory output'; end
    catch ME
        S.error = sprintf('%s: %s',ME.identifier,ME.message);
    end
end

% ========================================================================
function tOpts = makeTrajectoryOpts(opts,D,initialPose,off)
    tOpts = struct('camIDs',D.camIDs, 'initialPose',initialPose, ...
        'dataOverride',D, 'offsetsOverride',off, ...
        'nOuter',opts.trajectoryNOuter, ...
        'angleMidSigmaDeg',opts.trajectoryAngleMidSigmaDeg, ...
        'positionMidSigmaM',opts.trajectoryPositionMidSigmaM, ...
        'temporalHuberC',opts.trajectoryTemporalHuberC, ...
        'robustType',opts.robustType, 'robustC',opts.robustC, ...
        'saveOutput',false, 'verbose',false);
end

% ========================================================================
function A = anglesForTrajectory(pose,lureXYZ,Rz,Rx,Ry,nAng)
    nframes = size(pose,1);
    A = nan(nframes,nAng);
    for f = 1:nframes
        A(f,:) = mcMainAngles(pose(f,:),lureXYZ(f,:),Rz,Rx,Ry);
    end
end

% ========================================================================
function A = processAngles(Araw,dsp)
    A = Araw;
    if ~dsp.on, return; end
    for a = dsp.cols
        A(:,a) = mcFilterAngleSeries(A(:,a),isfinite(A(:,a)), ...
            dsp.K,dsp.nSigma,dsp.sgOrder,dsp.sgWin,ismember(a,dsp.unsignedCols));
    end
end

% ========================================================================
function rmse = poseFitRMSE(pose,D,off,Rz,Rx,Ry,sigmas,invPoses,intrin, ...
        ncams,nparts,partNames,robustOpt)
    nframes = size(pose,1);
    rmse = nan(nframes,1);
    for f = 1:nframes
        labF = squeeze(D.labels0(f,:,:,:));
        llhF = squeeze(D.llh(f,:,:));
        residual = mcMainResidual(pose(f,:),Rz,Rx,Ry,off,labF,llhF, ...
            sigmas,invPoses,intrin,ncams,nparts,partNames,robustOpt,D.camMaskStore(f,:));
        rmse(f) = sqrt(sum(residual.^2,'omitnan') / max(2*nnz(llhF>0.5),1));
    end
end

% ========================================================================
function S = aggregateAngles(samples,nominal)
    [~,nframes,nAng] = size(samples);
    devDeg = rad2deg(wrapToPi(samples - reshape(nominal,[1 nframes nAng])));
    medianDev = nan(nframes,nAng);
    loDev = nan(nframes,nAng);
    hiDev = nan(nframes,nAng);
    stdDev = nan(nframes,nAng);
    for f = 1:nframes
        for a = 1:nAng
            v = devDeg(:,f,a);
            v = v(isfinite(v));
            if isempty(v), continue; end
            medianDev(f,a) = median(v);
            q = prctile(v,[2.5 97.5]);
            loDev(f,a) = q(1);
            hiDev(f,a) = q(2);
            stdDev(f,a) = std(v);
        end
    end
    nominalDeg = rad2deg(nominal);
    S = struct('nominalDeg',nominalDeg, ...
        'medianDeg',nominalDeg+medianDev, ...
        'ciLoDeg',nominalDeg+loDev, ...
        'ciHiDeg',nominalDeg+hiDev, ...
        'ciHalfDeg',(hiDev-loDev)/2, ...
        'stdDeg',stdDev);
end

% ========================================================================
function off = offsetsFromLengths(L)
    % L = [eye_y nose_x nose_z ear_x ear_y ear_z neck_x neck_z spine_x spine_z]
    off = struct();
    off.l_eye     = [0;       L(1);  0];
    off.r_eye     = [0;      -L(1);  0];
    off.nose      = [L(2);    0;     L(3)];
    off.l_ear     = [L(4);    L(5);  L(6)];
    off.r_ear     = [L(4);   -L(5);  L(6)];
    off.neck_base = [L(7);    0;     L(8)];
    off.spine     = [L(9);    0;     L(10)];
end

% ========================================================================
function writeCsv(R, path)
    an = R.angleNames;
    T = table(R.absFrames,'VariableNames',{'absFrame'});
    for a = 1:numel(an)
        T.(sprintf('%s_nominal',an{a})) = R.nominalDeg(:,a);
        T.(sprintf('%s_median', an{a})) = R.medianDeg(:,a);
        T.(sprintf('%s_ciLo',   an{a})) = R.ciLoDeg(:,a);
        T.(sprintf('%s_ciHi',   an{a})) = R.ciHiDeg(:,a);
        T.(sprintf('%s_ciHalf', an{a})) = R.ciHalfDeg(:,a);
    end
    writetable(T, path);
end

% ========================================================================
function makeMainFigures(R, outDir, tag)
    an = R.angleNames; af = R.absFrames(:).';
    figure('Color','w','Position',[60 60 1500 820]);
    for a = 1:numel(an)
        subplot(2,3,a); hold on; grid on;
        lo = R.ciLoDeg(:,a).'; hi = R.ciHiDeg(:,a).';
        v = isfinite(lo)&isfinite(hi);
        if any(v), fill([af(v) fliplr(af(v))],[lo(v) fliplr(hi(v))],[0.3 0.5 0.9],'FaceAlpha',0.25,'EdgeColor','none'); end
        plot(af, R.medianDeg(:,a),'b-','LineWidth',1.4);
        plot(af, R.nominalDeg(:,a),'k--','LineWidth',0.8);
        xlabel('absolute frame'); ylabel([an{a} ' (deg)'],'Interpreter','none');
        title(sprintf('%s: CI half-width med=%.2f deg', an{a}, median(R.ciHalfDeg(:,a),'omitnan')),'Interpreter','none');
    end
    sgtitle(sprintf('%s: median (blue), nominal (dashed), 95%% CI band', R.datasetName),'Interpreter','none');
    saveas(gcf, fullfile(outDir, sprintf('mcMain_CI_bands%s.png', tag)));

    figure('Color','w','Position',[60 60 1100 480]); hold on; grid on;
    plot(af, R.ciHalfDeg,'LineWidth',1.3);
    legend(an,'Location','northeastoutside','Interpreter','none');
    xlabel('absolute frame'); ylabel('95% CI half-width (deg)');
    title(sprintf('%s: per-frame angle uncertainty', R.datasetName),'Interpreter','none');
    saveas(gcf, fullfile(outDir, sprintf('mcMain_CI_halfwidth%s.png', tag)));
end
