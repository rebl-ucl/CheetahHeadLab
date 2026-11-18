function R = reconTrajectoryBatch(datasetName, opts)
%RECONTRAJECTORYBATCH Robust bidirectional trajectory reconstruction.
%
% R = reconTrajectoryBatch(datasetName, opts)
%
% This is the bidirectional smoothing stage of CRBR, initialised by reconPerFrame.
% It starts from the robust per-frame CRBR poses, then repeatedly refits every
% interior frame to the original multi-camera labels with two terms:
%   1) the production Huber + cross-camera-consensus reprojection residual;
%   2) a robust constant-angular/linear-velocity trajectory prior formed from
%      BOTH neighbouring frames.
%
% Rotational priors are evaluated on SO(3), not on wrapped Euler yaw. A smooth
% 180-degree turn therefore has low cost; only changes in angular velocity are
% penalised. The temporal residual is Huberised so a genuine abrupt turn can
% override the prior when supported by the images.
%
% Outputs are deliberately separate from the existing reconstructions:
%   storedOptimRaw_trajectory.mat
%   trajectoryReconstruction/trajectory_diagnostics.csv
%   trajectoryReconstruction/trajectory_settings.mat
%
% Required inputs are the same files used by the Monte-Carlo loader, including
% storedOptimRaw.mat and triangulatedPointsRaw_multiview.mat. The latter
% supplies the production leave-one-out camera trust mask.
%
% opts fields (all optional):
%   dataRoot               <repo>/outputs
%   camIDs                 auto-detected
%   initialPose            [] (defaults to stored robust per-frame poses)
%   dataOverride           [] (validation-only mcLoadDataset-compatible struct)
%   excludeCamIDs          [] (validation-only whole-camera deletion)
%   offsetsOverride        [] (Monte-Carlo geometry override)
%   outSuffix              "_trajectory"
%   nOuter                 6
%   angleMidSigmaDeg       2.0  robust sigma for deviation from SO(3) midpoint
%   positionMidSigmaM      0.015
%   temporalHuberC         3.0
%   robustType/robustC     "huber" / 3 (production defaults)
%   dataWeight             1.0
%   maxIterations          250
%   maxFunctionEvaluations 5000
%   convergenceDeg         0.03 median pose change between outer iterations
%   saveOutput             true
%   verbose                true

if nargin < 2, opts = struct(); end

thisDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
mcDir = fullfile(repoRoot, 'CRBR', 'uncertainty');
addpath(thisDir, mcDir);

def = struct( ...
    'dataRoot', fullfile(repoRoot,'outputs'), ...
    'camIDs', [], 'outSuffix', "_trajectory", 'nOuter', 6, ...
    'initialPose', [], 'dataOverride', [], 'excludeCamIDs', [], ...
    'offsetsOverride', [], ...
    'angleMidSigmaDeg', 2.0, 'positionMidSigmaM', 0.015, ...
    'temporalHuberC', 3.0, 'robustType', 'huber', 'robustC', 3, ...
    'dataWeight', 1.0, 'maxIterations', 250, ...
    'maxFunctionEvaluations', 5000, 'convergenceDeg', 0.03, ...
    'saveOutput', true, 'verbose', true);
opts = fillDefaults(opts, def);

if isfolder(datasetName)
    parentFolder = char(datasetName);
    [~, datasetLabel] = fileparts(parentFolder);
else
    datasetLabel = char(datasetName);
    parentFolder = fullfile(opts.dataRoot, datasetLabel);
end
assert(isfolder(parentFolder), 'Dataset folder not found: %s', parentFolder);

if isempty(opts.camIDs) && isempty(opts.dataOverride)
    opts.camIDs = detectCameraIDs(parentFolder);
elseif isempty(opts.camIDs)
    opts.camIDs = opts.dataOverride.camIDs;
end
assert(~isempty(opts.camIDs), 'No cam*_labels_wide.csv files found in %s.', parentFolder);

if opts.verbose
    fprintf('\n=== Robust bidirectional trajectory reconstruction ===\n');
    fprintf('Dataset: %s\n', parentFolder);
    fprintf('Cameras: %s\n', mat2str(opts.camIDs));
end

if isempty(opts.dataOverride)
    D = mcLoadDataset(parentFolder, opts.camIDs);
else
    D = opts.dataOverride;
end
nFrames = D.nframes;
nCams = D.ncams;
nParts = D.nparts;
partNames = D.partNames;
sigmas = D.sigmas;
if isempty(opts.initialPose)
    nominal = D.xBarStore;
else
    nominal = opts.initialPose;
end
assert(size(nominal,1) == nFrames && size(nominal,2) == 12, ...
    'Expected an nFrames x 12 per-frame pose array.');
for cam = opts.excludeCamIDs(:).'
    j = find(D.camIDs == cam,1);
    if ~isempty(j), D.camMaskStore(:,j) = false; end
end

Rz = @(t)[cos(t) sin(t) 0; -sin(t) cos(t) 0; 0 0 1];
Rx = @(t)[1 0 0; 0 cos(t) sin(t); 0 -sin(t) cos(t)];
Ry = @(t)[cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
invPoses = cellfun(@(P) invert(P), D.poses, 'UniformOutput', false);
if isempty(opts.offsetsOverride)
    off = nominalOffsets(D.scaleHead0, D.scaleTorso0);
else
    off = opts.offsetsOverride;
end

robustOpt = struct('type', char(opts.robustType), 'c', opts.robustC, ...
    'consensus', struct('enable', ~strcmp(opts.robustType,'none'), ...
                        'C',4, 'floorSigma',2, 'minCams',3));
[LB, UB] = productionBounds(D.poses);
solveOpts = optimoptions(@lsqnonlin, 'Display','off', ...
    'StepTolerance',1e-5, 'OptimalityTolerance',1e-5, ...
    'MaxIterations',opts.maxIterations, ...
    'MaxFunctionEvaluations',opts.maxFunctionEvaluations);

pose = clampPose(nominal, LB, UB);
angleSigma = deg2rad(opts.angleMidSigmaDeg);
positionSigma = opts.positionMidSigmaM;
iteration = (1:opts.nOuter).';
medianChangeDeg = nan(opts.nOuter,1);
p95ChangeDeg = nan(opts.nOuter,1);
medianPositionChangeM = nan(opts.nOuter,1);
objectivePerFrame = nan(opts.nOuter,1);

for outer = 1:opts.nOuter
    oldPose = pose;
    newPose = oldPose;
    frameCost = nan(nFrames,1);

    % Endpoints remain data-determined. There is no symmetric acceleration
    % prior at an endpoint without inventing observations outside the clip.
    % Alternating forward/backward Gauss-Seidel updates converge faster than a
    % simultaneous Jacobi pass without introducing a permanent time direction.
    if mod(outer,2)==1
        frameOrder = 2:nFrames-1;
    else
        frameOrder = nFrames-1:-1:2;
    end
    for f = frameOrder
        target = midpointTarget(newPose(f-1,:), newPose(f+1,:), Rz,Rx,Ry);
        labF = squeeze(D.labels0(f,:,:,:));
        llhF = squeeze(D.llh(f,:,:));
        camMask = D.camMaskStore(f,:);
        dataFun = @(p) mcMainResidual(p, Rz,Rx,Ry, off, labF, llhF, ...
            sigmas, invPoses, D.intrin, nCams, nParts, partNames, robustOpt, camMask);
        obj = @(p) trajectoryResidual(p, dataFun, target, Rz,Rx,Ry, ...
            angleSigma, positionSigma, opts.temporalHuberC, opts.dataWeight);

        starts = {clampPose(oldPose(f,:),LB,UB), clampPose(nominal(f,:),LB,UB)};
        bestP = starts{1}; bestCost = inf;
        for s = 1:numel(starts)
            try
                [candidate, cost] = lsqnonlin(obj, starts{s}, LB, UB, solveOpts);
                if isfinite(cost) && cost < bestCost
                    bestP = candidate;
                    bestCost = cost;
                end
            catch ME
                if opts.verbose
                    warning('reconTrajectoryBatch:frameSolve', ...
                        'Frame %d, start %d failed: %s', f, s, ME.message);
                end
            end
        end
        newPose(f,:) = bestP;
        frameCost(f) = bestCost;
    end

    pose = newPose;
    [rotChange, posChange] = poseChange(oldPose, pose, Rz,Rx,Ry);
    medianChangeDeg(outer) = median(rotChange,'omitnan');
    p95ChangeDeg(outer) = prctile(rotChange,95);
    medianPositionChangeM(outer) = median(posChange,'omitnan');
    objectivePerFrame(outer) = median(frameCost,'omitnan');
    if opts.verbose
        fprintf(['Outer %d/%d: median rotation update %.3f deg; p95 %.3f deg; ' ...
            'median position update %.4f m\n'], outer, opts.nOuter, ...
            medianChangeDeg(outer), p95ChangeDeg(outer), medianPositionChangeM(outer));
    end
    if medianChangeDeg(outer) < opts.convergenceDeg
        iteration = iteration(1:outer);
        medianChangeDeg = medianChangeDeg(1:outer);
        p95ChangeDeg = p95ChangeDeg(1:outer);
        medianPositionChangeM = medianPositionChangeM(1:outer);
        objectivePerFrame = objectivePerFrame(1:outer);
        break
    end
end

diagnostics = table(iteration, medianChangeDeg, p95ChangeDeg, ...
    medianPositionChangeM, objectivePerFrame, ...
    'VariableNames', {'iteration','median_rotation_update_deg', ...
    'p95_rotation_update_deg','median_position_update_m','median_objective_per_frame'});

xBarRawStore = pose;
xBarStore = pose;
absFrames = D.absFrames(:);
seedSourceLog = repmat("robust_bidirectional", nFrames, 1);
visible = all(isfinite(D.labels0),4) & D.llh > 0.5;
visible = visible & reshape(D.camMaskStore,[nFrames 1 nCams]);
nValidMVperFrame = sum(sum(visible,3)>=2,2);
reconMeta = struct('reconMode',"robust_bidirectional", 'source',"raw", ...
    'temporal',true, 'bidirectional',true, 'rotationManifold',"SO3", ...
    'trajectoryPrior',"robust_constant_velocity", ...
    'robustType',string(opts.robustType), 'robustC',opts.robustC, ...
    'angleMidSigmaDeg',opts.angleMidSigmaDeg, ...
    'positionMidSigmaM',opts.positionMidSigmaM, ...
    'temporalHuberC',opts.temporalHuberC, 'nOuter',height(diagnostics), ...
    'absFrameStart',absFrames(1), 'absFrameEnd',absFrames(end), ...
    'nFrames',nFrames, 'dataset',string(datasetLabel));

outDir = fullfile(parentFolder, 'trajectoryReconstruction');
if opts.saveOutput
    if ~isfolder(outDir), mkdir(outDir); end
    outFile = fullfile(parentFolder, "storedOptimRaw" + string(opts.outSuffix) + ".mat");
    save(outFile, 'xBarRawStore','xBarStore','seedSourceLog', ...
        'nValidMVperFrame','absFrames','reconMeta');
    writetable(diagnostics, fullfile(outDir,'trajectory_diagnostics.csv'));
    settings = opts;
    if isfield(settings,'dataOverride'), settings.dataOverride = []; end
    if isfield(settings,'offsetsOverride') && ~isempty(settings.offsetsOverride)
        settings.offsetsOverride = "provided_at_runtime";
    end
    if isfield(settings,'initialPose') && ~isempty(settings.initialPose)
        settings.initialPoseSize = size(settings.initialPose);
        settings.initialPose = [];
    end
    save(fullfile(outDir,'trajectory_settings.mat'),'settings','reconMeta');
    if opts.verbose, fprintf('Saved %s\n', outFile); end
end

R = struct('dataset',string(datasetLabel), 'parentFolder',string(parentFolder), ...
    'pose',pose, 'nominalPose',nominal, 'absFrames',absFrames, ...
    'diagnostics',diagnostics, 'reconMeta',reconMeta, 'outDir',string(outDir));
end

function r = trajectoryResidual(p, dataFun, target, Rz,Rx,Ry, angleSigma, posSigma, huberC, dataWeight)
dataR = sqrt(dataWeight) * dataFun(p);
[Rh,Rb] = rotationsFromPose(p,Rz,Rx,Ry);
headRot = so3Log(target.Rh' * Rh) / angleSigma;
bodyRot = so3Log(target.Rb' * Rb) / angleSigma;
headPos = (p(4:6).' - target.headPos) / posSigma;
bodyPos = (p(10:12).' - target.bodyPos) / posSigma;
temp = [headRot; bodyRot; headPos; bodyPos];
temp = huberPseudoResidual(temp, huberC);
r = [dataR(:); temp(:)];
end

function target = midpointTarget(a,b,Rz,Rx,Ry)
[Rah,Rab] = rotationsFromPose(a,Rz,Rx,Ry);
[Rbh,Rbb] = rotationsFromPose(b,Rz,Rx,Ry);
target.Rh = rotationMidpoint(Rah,Rbh);
target.Rb = rotationMidpoint(Rab,Rbb);
target.headPos = 0.5*(a(4:6).' + b(4:6).');
target.bodyPos = 0.5*(a(10:12).' + b(10:12).');
end

function [Rh,Rb] = rotationsFromPose(p,Rz,Rx,Ry)
Rh = Rz(p(3))*Rx(p(1))*Ry(p(2));
Rb = Rz(p(9))*Rx(p(7))*Ry(p(8));
end

function Rm = rotationMidpoint(Ra,Rb)
v = so3Log(Ra' * Rb);
Rm = Ra * so3Exp(0.5*v);
end

function v = so3Log(R)
R = projectRotation(R);
c = max(-1,min(1,(trace(R)-1)/2));
theta = acos(c);
if theta < 1e-8
    v = 0.5 * [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
elseif pi-theta < 1e-5
    [V,D] = eig(R);
    [~,k] = min(abs(diag(D)-1));
    axis = real(V(:,k));
    axis = axis/max(norm(axis),eps);
    skewHint = [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
    if dot(axis,skewHint) < 0, axis = -axis; end
    v = theta*axis;
else
    v = theta/(2*sin(theta)) * ...
        [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
end
end

function R = so3Exp(v)
theta = norm(v);
K = skew(v);
if theta < 1e-8
    R = eye(3) + K + 0.5*K*K;
else
    R = eye(3) + (sin(theta)/theta)*K + ((1-cos(theta))/theta^2)*(K*K);
end
R = projectRotation(R);
end

function K = skew(v)
K = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

function R = projectRotation(A)
[U,~,V] = svd(A);
R = U*diag([1 1 det(U*V')])*V';
end

function y = huberPseudoResidual(x,c)
a = abs(x);
y = x;
m = a > c;
y(m) = sign(x(m)).*sqrt(max(2*c*a(m)-c^2,0));
end

function [rotDeg,posM] = poseChange(a,b,Rz,Rx,Ry)
n = size(a,1);
rotDeg = nan(2*n,1);
posM = nan(2*n,1);
for i=1:n
    [A1,A2] = rotationsFromPose(a(i,:),Rz,Rx,Ry);
    [B1,B2] = rotationsFromPose(b(i,:),Rz,Rx,Ry);
    rotDeg(2*i-1:2*i) = rad2deg([norm(so3Log(A1'*B1)); norm(so3Log(A2'*B2))]);
    posM(2*i-1:2*i) = [norm(a(i,4:6)-b(i,4:6)); norm(a(i,10:12)-b(i,10:12))];
end
end

function off = nominalOffsets(scaleHead,scaleTorso)
off = struct();
off.l_eye = scaleHead*[0;0.03;0];
off.r_eye = scaleHead*[0;-0.03;0];
off.nose = scaleHead*[0.055;0;-0.055];
off.l_ear = scaleHead*[-0.0725;0.0608;0.0441];
off.r_ear = scaleHead*[-0.0725;-0.0608;0.0441];
off.neck_base = scaleTorso*[0.04;0;0.10];
off.spine = scaleTorso*[-0.33;0;0.10];
end

function [LB,UB] = productionBounds(poses)
xc = cellfun(@(P) P.Translation(1),poses);
yc = cellfun(@(P) P.Translation(2),poses);
zc = cellfun(@(P) P.Translation(3),poses);
LB = [-pi/5,-pi/3,-pi,min(xc)-5,min(yc)-2,min(zc)-0.5, ...
      -pi/5,-pi/3,-pi,min(xc)-5,min(yc)-2,min(zc)-0.5];
UB = [ pi/5, pi/3, pi,max(xc)+5,max(yc)+2,max(zc)+1.5, ...
       pi/5, pi/3, pi,max(xc)+5,max(yc)+2,max(zc)+1.5];
end

function p = clampPose(p,LB,UB)
p = min(max(p,LB),UB);
end

function ids = detectCameraIDs(folder)
f = dir(fullfile(folder,'cam*_labels_wide.csv'));
ids = [];
for i=1:numel(f)
    tok = regexp(f(i).name,'^cam(\d+)_labels_wide\.csv$','tokens','once');
    if ~isempty(tok), ids(end+1) = str2double(tok{1}); end %#ok<AGROW>
end
ids = unique(ids,'sorted');
end

function opts = fillDefaults(opts,def)
f = fieldnames(def);
for i=1:numel(f)
    if ~isfield(opts,f{i}) || isempty(opts.(f{i})), opts.(f{i})=def.(f{i}); end
end
end
