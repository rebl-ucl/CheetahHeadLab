function R = yawSpatialLookAheadAnalysis(opts)
% yawSpatialLookAheadAnalysis  Eye-fixed spatial look-ahead on the lure path.
%
% R = yawSpatialLookAheadAnalysis(opts)
%
% Restricted to the five steady-lateral pursuit runs. For each valid frame t:
%   1. cast the current horizontal head-forward ray from the current eye;
%   2. intersect it continuously with the piecewise-linear full lure path;
%   3. choose the first forward intersection P_H(t); and
%   4. calculate d(t) = s(P_H(t)) - s(L(t)), where s is cumulative 3-D
%      arc length along the time-oriented lure trajectory.
%
% Positive d is ahead of the current lure along its path; negative d is
% behind. This function writes numerical tables only. The manuscript
% time series are exported by yawSpatialLookAheadMonteCarloAnalysis.

arguments
    opts struct = struct()
end

thisDir = fileparts(mfilename('fullpath'));
pursuitDir = fileparts(thisDir);
repoRoot = fileparts(thisDir);
opts = applyDefaults(opts, struct( ...
    'dataRoot', fullfile(repoRoot, 'outputs', 'analysis_inputs', 'yaw_tracking'), ...
    'outputDir', fullfile(repoRoot, 'outputs', 'yaw', 'spatial'), ...
    'displaySmoothSec', 0.075, ...
    'extendEndpoint', true, 'endpointWindowSec', 0.100, ...
    'zeroToleranceM', 0.05, ...
    'showFigures', false, 'saveFigures', false, ...
    'figureResolution', 600));

assert(~opts.saveFigures,'release:nonManuscriptPlot', ...
    ['Nominal-only spatial plots are not in the manuscript. Use ' ...
     'yawSpatialLookAheadMonteCarloAnalysis for Figure S15b,c.']);

assert(isscalar(opts.displaySmoothSec) && isfinite(opts.displaySmoothSec) && ...
    opts.displaySmoothSec >= 0, 'displaySmoothSec must be non-negative.');
assert(isscalar(opts.zeroToleranceM) && isfinite(opts.zeroToleranceM) && ...
    opts.zeroToleranceM >= 0, 'zeroToleranceM must be non-negative.');

if ~exist(opts.outputDir, 'dir')
    mkdir(opts.outputDir);
end

C = steadyLateralConfiguration();
fprintf('\n==== yawSpatialLookAheadAnalysis ====\n');
fprintf('Steady-lateral datasets: %d\n', height(C));
fprintf('Output directory: %s\n', opts.outputDir);

runs = cell(height(C),1);
frameRows = cell(height(C),1);
summaryRows = cell(height(C),1);
for i = 1:height(C)
    dataset = char(C.dataset(i));
    datasetFolder = fullfile(opts.dataRoot, dataset);
    resultPath = fullfile(datasetFolder, 'lureAnalysis_results_trajectory.mat');
    assert(isfile(resultPath), 'Missing trajectory result: %s', resultPath);
    S = load(resultPath, 'runResults');
    assert(isfield(S,'runResults'), '%s lacks runResults.', resultPath);

    Q = estimateRun(S.runResults, datasetFolder, C.dataset(i), C.animal(i), opts);
    runs{i} = Q;
    frameRows{i} = Q.frameTable;
    summaryRows{i} = Q.summaryTable;
    fprintf(['  %-30s intersections=%3d/%3d, median d=%7.2f m, ' ...
        'rho(d,range)=%6.2f\n'], dataset, Q.summaryTable.nIntersectionFrames, ...
        Q.summaryTable.nAnalysisFrames, Q.summaryTable.medianLeadDistanceM, ...
        Q.summaryTable.spearmanLeadVsHeadLureDistance);
end

frameTable = vertcat(frameRows{:});
summaryTable = vertcat(summaryRows{:});
writetable(frameTable, fullfile(opts.outputDir, ...
    'yaw_spatial_lookahead_by_frame.csv'));
writetable(summaryTable, fullfile(opts.outputDir, ...
    'yaw_spatial_lookahead_summary_by_run.csv'));

% Nominal-only time-series figures are not used in the manuscript.
writeReadme(C, opts);

R = struct('configuration',C, 'runs',{runs}, 'frameTable',frameTable, ...
    'summaryTable',summaryTable, 'opts',opts);
save(fullfile(opts.outputDir, 'yaw_spatial_lookahead_analysis.mat'), 'R');
fprintf('Saved steady-lateral spatial look-ahead analysis.\n');
end

function C = steadyLateralConfiguration()
dataset = ["julesFlick1_20190309"; "julesTopRun11_20170829"; ...
    "julesTopRun1_20171209"; "lilyRun_20190309"; "phantomRun_20190303"];
animal = ["Jules"; "Jules"; "Jules"; "Lily"; "Phantom"];
C = table(dataset,animal);
end

function Q = estimateRun(A, datasetFolder, dataset, animal, opts)
required = {'commonFrames','theta_h_yaw','lambda_h_yaw','r'};
for i = 1:numel(required)
    assert(isfield(A,required{i}), 'runResults lacks %s for %s.', ...
        required{i}, dataset);
end

absFrame = double(A.commonFrames(:));
thetaHead = double(A.theta_h_yaw(:));
lambdaYaw = double(A.lambda_h_yaw(:));
headLureDistance = double(A.r(:));
N = numel(absFrame);
assert(all([numel(thetaHead),numel(lambdaYaw),numel(headLureDistance)] == N), ...
    'Yaw/range vectors do not match commonFrames for %s.', dataset);

if isfield(A,'fs') && ~isempty(A.fs) && isfinite(A.fs(1)) && A.fs(1)>0
    fs = double(A.fs(1));
else
    fs = 90;
    warning('yawSpatialLookAhead:fpsFallback', ...
        '%s lacks a valid runResults.fs; using 90 Hz.', dataset);
end
if isfield(A,'lambda_h_pitch') && numel(A.lambda_h_pitch)==N
    lambdaPitch = double(A.lambda_h_pitch(:));
else
    lambdaPitch = zeros(N,1);
end
if isfield(A,'yawAnalysisMask') && numel(A.yawAnalysisMask)==N
    analysisMask = logical(A.yawAnalysisMask(:));
else
    analysisMask = true(N,1);
end
if isfield(A,'isInterp') && numel(A.isInterp)==N
    poseIsInterpolated = logical(A.isInterp(:));
else
    poseIsInterpolated = false(N,1);
end

[pathFrame,pathXYZ,pathIsInterpolated] = loadFullLureTrajectory(datasetFolder);
segmentLength = sqrt(sum(diff(pathXYZ,1,1).^2,2));
pathS = [0; cumsum(segmentLength)];

% Recover current eye position by inverting the stored eye-to-current-lure
% spherical vector. Full and cropped lure coordinates agree on these frames.
[hasCurrentLure,currentLoc] = ismember(absFrame,pathFrame);
currentLureXYZ = nan(N,3);
currentLureS = nan(N,1);
currentLureXYZ(hasCurrentLure,:) = pathXYZ(currentLoc(hasCurrentLure),:);
currentLureS(hasCurrentLure) = pathS(currentLoc(hasCurrentLure));
[vx,vy,vz] = sph2cart(lambdaYaw,lambdaPitch,headLureDistance);
eyeXYZ = currentLureXYZ-[vx,vy,vz];

currentValid = analysisMask & all(isfinite(eyeXYZ),2) & ...
    isfinite(thetaHead) & isfinite(currentLureS) & isfinite(headLureDistance);
leadDistance = nan(N,1);
intersectionXYZ = nan(N,3);
intersectionS = nan(N,1);
rayDistance = nan(N,1);
intersectionCount = zeros(N,1);
multipleIntersections = false(N,1);
chosenSegmentStartFrame = nan(N,1);
chosenSegmentEndFrame = nan(N,1);
intersectionUsesInterpolatedLure = false(N,1);
intersectionUsesEndpointExtension = false(N,1);
terminalDirection=spatialTerminalDirection(pathFrame,pathXYZ,fs,opts.endpointWindowSec);

for k = find(currentValid).'
    origin = eyeXYZ(k,1:2);
    direction = [cos(thetaHead(k)),sin(thetaHead(k))];
    H = rayPathIntersections(origin,direction,pathFrame,pathXYZ,pathS, ...
        pathIsInterpolated);
    intersectionCount(k) = H.n;
    multipleIntersections(k) = H.n>1;
    if H.n==0
        if opts.extendEndpoint
            [ld,rho,point]=spatialForwardExtension(origin,direction,pathXYZ(end,:), ...
                terminalDirection,pathS(end),currentLureS(k));
            if isfinite(ld)
                leadDistance(k)=ld; intersectionXYZ(k,:)=point;
                intersectionS(k)=currentLureS(k)+ld; rayDistance(k)=rho;
                intersectionCount(k)=1; intersectionUsesEndpointExtension(k)=true;
                chosenSegmentStartFrame(k)=pathFrame(end);
            end
        end
        continue
    end
    minRay = min(H.rayDistanceM);
    tied = find(abs(H.rayDistanceM-minRay) <= 1e-8);
    if numel(tied)>1
        [~,j] = min(abs(H.pathS(tied)-currentLureS(k)));
        chosen = tied(j);
    else
        chosen = tied;
    end
    leadDistance(k) = H.pathS(chosen)-currentLureS(k);
    intersectionXYZ(k,:) = H.xyz(chosen,:);
    intersectionS(k) = H.pathS(chosen);
    rayDistance(k) = H.rayDistanceM(chosen);
    chosenSegmentStartFrame(k) = H.segmentStartFrame(chosen);
    chosenSegmentEndFrame(k) = H.segmentEndFrame(chosen);
    intersectionUsesInterpolatedLure(k) = H.usesInterpolatedLure(chosen);
end

smoothFrames = max(1,round(opts.displaySmoothSec*fs));
if mod(smoothFrames,2)==0, smoothFrames=smoothFrames+1; end
if opts.displaySmoothSec>0
    leadSmooth = movmedian(leadDistance,smoothFrames,'omitmissing', ...
        'Endpoints','shrink');
else
    leadSmooth = leadDistance;
end
leadSmooth(~isfinite(leadDistance)) = NaN;
timeSec = (absFrame-absFrame(1))/fs;

frameTable = table(repmat(dataset,N,1),repmat(animal,N,1),absFrame,timeSec, ...
    repmat(fs,N,1),analysisMask,poseIsInterpolated,eyeXYZ(:,1),eyeXYZ(:,2), ...
    eyeXYZ(:,3),headLureDistance,currentLureS,leadDistance,leadSmooth, ...
    intersectionXYZ(:,1),intersectionXYZ(:,2),intersectionXYZ(:,3), ...
    intersectionS,rayDistance,intersectionCount,multipleIntersections, ...
    chosenSegmentStartFrame,chosenSegmentEndFrame, ...
    intersectionUsesInterpolatedLure, ...
    'VariableNames',{'dataset','animal','absFrame','timeSec','fps', ...
    'analysisMask','isInterpolatedPose','eyeX_m','eyeY_m','eyeZ_m', ...
    'headOriginToLureDistanceM','currentLurePathS_m','leadDistanceM', ...
    'leadDistanceDisplaySmoothM','intersectionX_m','intersectionY_m', ...
    'intersectionZ_m','intersectionPathS_m','horizontalRayDistanceM', ...
    'intersectionCount','multipleIntersections','segmentStartAbsFrame', ...
    'segmentEndAbsFrame','intersectionUsesInterpolatedLure'});

frameTable.intersectionUsesEndpointExtension=intersectionUsesEndpointExtension;
summaryTable = summariseRun(frameTable,dataset,animal,opts.zeroToleranceM);
summaryTable.nEndpointExtensionFrames=nnz(intersectionUsesEndpointExtension);
summaryTable.nFullLureFrames = numel(pathFrame);
summaryTable.fullLurePathLengthM = pathS(end);
Q = struct('dataset',dataset,'animal',animal,'fs',fs, ...
    'frameTable',frameTable,'summaryTable',summaryTable);
end

function H = rayPathIntersections(origin,direction,frame,xyz,pathS,isInterp)
% Intersections of a forward horizontal ray with a time-ordered 3-D path.
% The crossing is solved in XY; pathS and intersection Z interpolate along
% the same 3-D segment.
H = struct('n',0,'rayDistanceM',zeros(0,1),'pathS',zeros(0,1), ...
    'xyz',zeros(0,3),'segmentStartFrame',zeros(0,1), ...
    'segmentEndFrame',zeros(0,1),'usesInterpolatedLure',false(0,1));
dn = norm(direction);
if dn<1e-12 || size(xyz,1)<2, return; end
d = direction/dn;

for k = 1:size(xyz,1)-1
    A = xyz(k,1:2); B = xyz(k+1,1:2); e = B-A;
    if any(~isfinite([A,B,frame(k),frame(k+1),pathS(k),pathS(k+1)]))
        continue
    end
    q = A-origin;
    denominator = cross2(d,e);
    scaleTol = 1e-10*max([1,norm(e),norm(q)]);
    if abs(denominator)>scaleTol
        rho = cross2(q,e)/denominator;
        alpha = cross2(q,d)/denominator;
        if rho < -scaleTol || alpha < -scaleTol || alpha > 1+scaleTol
            continue
        end
        rho = max(rho,0); alpha = min(max(alpha,0),1);
    else
        % Parallel non-collinear segments do not intersect. For a collinear
        % segment, use its nearest forward point.
        if abs(cross2(q,d))>scaleTol, continue; end
        rhoA = dot(A-origin,d); rhoB = dot(B-origin,d);
        if max(rhoA,rhoB)<-scaleTol, continue; end
        if min(rhoA,rhoB)>=0
            if rhoA<=rhoB, rho=rhoA; alpha=0; else, rho=rhoB; alpha=1; end
        else
            rho = 0;
            alpha = -rhoA/(rhoB-rhoA);
            alpha = min(max(alpha,0),1);
        end
    end

    hitXYZ = xyz(k,:)+alpha*(xyz(k+1,:)-xyz(k,:));
    hitS = pathS(k)+alpha*(pathS(k+1)-pathS(k));
    % Adjacent segments duplicate a shared vertex. Preserve spatially
    % coincident hits only when they represent different path branches.
    duplicate = false;
    for j = 1:H.n
        if norm(H.xyz(j,:)-hitXYZ)<1e-8 && abs(H.pathS(j)-hitS)<1e-8
            duplicate = true; break
        end
    end
    if duplicate, continue; end
    H.n = H.n+1;
    H.rayDistanceM(end+1,1) = rho;
    H.pathS(end+1,1) = hitS;
    H.xyz(end+1,:) = hitXYZ;
    H.segmentStartFrame(end+1,1) = frame(k);
    H.segmentEndFrame(end+1,1) = frame(k+1);
    H.usesInterpolatedLure(end+1,1) = isInterp(k) || isInterp(k+1);
end

if H.n>1
    [~,order] = sort(H.rayDistanceM);
    fields = {'rayDistanceM','pathS','xyz','segmentStartFrame', ...
        'segmentEndFrame','usesInterpolatedLure'};
    for i = 1:numel(fields)
        H.(fields{i}) = H.(fields{i})(order,:);
    end
end
end

function z = cross2(a,b)
z = a(1)*b(2)-a(2)*b(1);
end

function T = summariseRun(F,dataset,animal,zeroToleranceM)
valid = isfinite(F.leadDistanceM) & isfinite(F.headOriginToLureDistanceM);
d = F.leadDistanceM(valid);
r = F.headOriginToLureDistanceM(valid);
t = F.timeSec(valid);
if isempty(d)
    med=NaN; q25=NaN; q75=NaN; fracAhead=NaN; fracAligned=NaN;
    fracBehind=NaN; rho=NaN; slopeRange=NaN; slopeTime=NaN;
    firstMed=NaN; lastMed=NaN; firstRange=NaN; lastRange=NaN;
else
    med=median(d); q25=percentile(d,25); q75=percentile(d,75);
    fracAhead=mean(d>zeroToleranceM);
    fracAligned=mean(abs(d)<=zeroToleranceM);
    fracBehind=mean(d<-zeroToleranceM);
    rho=spearmanCorrelation(d,r);
    slopeRange=linearSlope(r,d);
    slopeTime=linearSlope(t,d);
    n=numel(d); third=max(1,floor(n/3));
    firstMed=median(d(1:third)); lastMed=median(d(end-third+1:end));
    firstRange=median(r(1:third)); lastRange=median(r(end-third+1:end));
end
T = table(dataset,animal,nnz(F.analysisMask),nnz(valid), ...
    nnz(~isfinite(F.leadDistanceM) & F.analysisMask)/max(1,nnz(F.analysisMask)), ...
    mean(F.multipleIntersections(valid)), ...
    med,q25,q75,q75-q25,fracAhead,fracAligned,fracBehind, ...
    median(r,'omitnan'),rho,slopeRange,slopeTime,firstMed,lastMed, ...
    lastMed-firstMed,firstRange,lastRange,lastRange-firstRange, ...
    'VariableNames',{'dataset','animal','nAnalysisFrames','nIntersectionFrames', ...
    'fractionNoIntersection','fractionMultipleIntersections', ...
    'medianLeadDistanceM','q25LeadDistanceM','q75LeadDistanceM', ...
    'iqrLeadDistanceM','fractionAhead','fractionCurrentLureAligned', ...
    'fractionBehind','medianHeadLureDistanceM', ...
    'spearmanLeadVsHeadLureDistance','linearLeadPerDistanceSlope', ...
    'linearLeadTimeSlopeMPerSec','firstThirdMedianLeadM', ...
    'lastThirdMedianLeadM','lastMinusFirstThirdLeadM', ...
    'firstThirdMedianHeadLureDistanceM','lastThirdMedianHeadLureDistanceM', ...
    'lastMinusFirstThirdHeadLureDistanceM'});
end

function [frame,xyz,isInterp] = loadFullLureTrajectory(datasetFolder)
path=fullfile(datasetFolder,'3DLurePos_full.mat');
assert(isfile(path),'Missing full lure trajectory: %s',path);
S=load(path,'Lure3D');
assert(isfield(S,'Lure3D') && isfield(S.Lure3D,'Index') && ...
    isfield(S.Lure3D,'Mark'),'Expected Lure3D.Index/Mark in %s.',path);
L=S.Lure3D; frame=double(L.Index(:)); xyz=nan(numel(frame),3);
for i=1:numel(frame)
    mark=L.Mark{i};
    if ~isempty(mark) && numel(mark)>=3, xyz(i,:)=double(mark(1:3)); end
end
if isfield(L,'IsInterpolated') && numel(L.IsInterpolated)==numel(frame)
    isInterp=logical(L.IsInterpolated(:));
else
    isInterp=false(numel(frame),1);
end
valid=isfinite(frame)&all(isfinite(xyz),2);
frame=frame(valid); xyz=xyz(valid,:); isInterp=isInterp(valid);
[frame,order]=sort(frame); xyz=xyz(order,:); isInterp=isInterp(order);
assert(~isempty(frame),'No finite lure positions in %s.',path);
assert(numel(unique(frame))==numel(frame),'Duplicate lure frames in %s.',path);
end

function rho = spearmanCorrelation(x,y)
valid=isfinite(x)&isfinite(y); x=x(valid); y=y(valid);
if numel(x)<3, rho=NaN; return; end
rx=averageRanks(x); ry=averageRanks(y);
c=corrcoef(rx,ry);
if numel(c)<4, rho=NaN; else, rho=c(1,2); end
end

function r = averageRanks(x)
[sorted,order]=sort(x(:)); r=zeros(size(x(:))); n=numel(x); i=1;
while i<=n
    j=i;
    while j<n && sorted(j+1)==sorted(i), j=j+1; end
    r(order(i:j))=(i+j)/2;
    i=j+1;
end
end

function slope = linearSlope(x,y)
valid=isfinite(x)&isfinite(y); x=x(valid); y=y(valid);
if numel(x)<2 || range(x)==0, slope=NaN; return; end
p=polyfit(x,y,1); slope=p(1);
end

function p = percentile(x,q)
x=sort(x(isfinite(x)));
if isempty(x), p=NaN; return; end
u=1+(numel(x)-1)*q/100; lo=floor(u); hi=ceil(u);
if lo==hi, p=x(lo); else, p=x(lo)+(u-lo)*(x(hi)-x(lo)); end
end

function writeReadme(C,opts)
path=fullfile(opts.outputDir,'README.txt');
fid=fopen(path,'w'); assert(fid>=0,'Could not write %s.',path);
fprintf(fid,'Endpoint continuation enabled: %d; trailing direction window %.3f s.\n',opts.extendEndpoint,opts.endpointWindowSec);
fprintf(fid,'When enabled, missing observed-path intersections use a forward horizontal ray from the final target position at constant height.\n');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'STEADY-LATERAL SPATIAL LOOK-AHEAD\n\n');
fprintf(fid,['For each frame, the current horizontal head-forward ray is ' ...
    'intersected continuously with the full lure polyline.\n']);
fprintf(fid,['The first forward intersection is selected. Signed d is its ' ...
    'cumulative 3-D path arc length minus the current lure arc length.\n']);
fprintf(fid,['Positive d is ahead; zero is current-lure alignment; negative d ' ...
    'is behind along the time-oriented path.\n\n']);
fprintf(fid,['The %.0f ms moving median is a visual guide only. The right axis ' ...
    'is current 3-D head-origin-to-lure distance from runResults.r.\n\n'], ...
    1000*opts.displaySmoothSec);
fprintf(fid,['Multiple ray/path crossings are reported. The nearest forward ' ...
    'intersection is primary; coincident branches are tie-broken to the ' ...
    'current path location.\n\n']);
fprintf(fid,'Included datasets (%d):\n',height(C));
for i=1:height(C), fprintf(fid,'  %s\n',C.dataset(i)); end
end

function opts = applyDefaults(opts,defaults)
names=fieldnames(defaults);
for i=1:numel(names)
    if ~isfield(opts,names{i}) || isempty(opts.(names{i}))
        opts.(names{i})=defaults.(names{i});
    end
end
end
