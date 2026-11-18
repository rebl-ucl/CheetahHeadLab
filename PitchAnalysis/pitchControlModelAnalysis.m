function A = pitchControlModelAnalysis(datasetList, opts)
% pitchControlModelAnalysis  Direct tests of the functional pitch controller.
%
% The measured sagittal-plane variables are
%   lambda = head-to-lure elevation (reference)
%   thetaH = world-frame head elevation (controlled output)
%   thetaT = world-frame torso pitch (disturbance)
%   q      = thetaH-thetaT (head-on-torso pitch)
%   alpha  = thetaH-lambda (tracking error)
%
% Three pre-specified finite-change models are fitted separately to every
% pursuit sequence and then averaged run-equal within animal:
%   disturbance rejection: dThetaH = bLambda*dLambda + bTorso*dThetaT
%   counter-rotation:       dQ      = -gD*dThetaT + gLambda*dLambda
%   delayed correction:    dQ(t+tau) = -kAlpha*alpha(t)
% The adjusted delayed model also controls the torso and reference changes
% occurring during the future response window. Positive kAlpha means that
% subsequent head-on-torso movement tends to reduce the current pitch error.

arguments
    datasetList
    opts struct = struct()
end

thisDir = fileparts(mfilename('fullpath'));
lureDir = thisDir;
repoRoot = fileparts(thisDir);
opts = applyDefaults(opts, struct( ...
    'dataRoot', fullfile(repoRoot,'outputs','analysis_inputs'), ...
    'outputDir', fullfile(repoRoot,'outputs','pitch','mechanism'), ...
    'resultFileName', 'lureAnalysis_results_trajectory.mat', ...
    'frameRanges', {cell(0,2)}, ...
    'frameRangeMode', 'absFrame', ...
    'excludeInterpolated', false, ...
    'smoothSec', 0.05, ...
    'changeWindowSec', 0.05, ...
    'primaryCorrectionDelaySec', 0.05, ...
    'correctionDelayGridSec', (0.01:0.01:0.15), ...
    'minPairsPerRun', 20, ...
    'representativeDataset', 'pitch_tracking/phantomTopRun1_20170903', ...
    'representativeSmoothSec', 0, ...
    'saveFigures', true, ...
    'showFigures', false, ...
    'figureResolution', 600));

opts.dataRoot = char(opts.dataRoot);
opts.outputDir = char(opts.outputDir);
opts.resultFileName = char(opts.resultFileName);
opts.frameRangeMode = validatestring(opts.frameRangeMode, {'absFrame','alignedRow'});
assert(opts.smoothSec >= 0, 'smoothSec cannot be negative.');
assert(opts.changeWindowSec > 0, 'changeWindowSec must be positive.');
assert(opts.primaryCorrectionDelaySec >= 0, ...
    'primaryCorrectionDelaySec cannot be negative.');
assert(all(opts.correctionDelayGridSec >= 0), ...
    'correctionDelayGridSec cannot contain negative delays.');
assert(opts.minPairsPerRun >= 4, 'minPairsPerRun must be at least four.');
if ~exist(opts.outputDir,'dir'), mkdir(opts.outputDir); end

datasetList = normaliseList(datasetList);
assert(~isempty(datasetList), 'datasetList cannot be empty.');
requestedDelays = unique([opts.primaryCorrectionDelaySec, ...
    opts.correctionDelayGridSec(:).']);

runRows = {};
fprintf('\n==== pitchControlModelAnalysis ====\n');
fprintf('Finite-change window %.0f ms; primary correction delay %.0f ms.\n', ...
    1000*opts.changeWindowSec, 1000*opts.primaryCorrectionDelaySec);

for ri = 1:numel(datasetList)
    key = normaliseKey(datasetList{ri});
    [group,dataset] = splitDatasetKey(key);
    animal = inferAnimal(dataset);
    resultPath = fullfile(opts.dataRoot,group,dataset,opts.resultFileName);
    assert(isfile(resultPath), 'Missing result for %s: %s',key,resultPath);
    S = load(resultPath,'runResults'); R = S.runResults;
    required = {'theta_h_pitch','theta_t_pitch','lambda_h_pitch', ...
        'pitchAnalysisMask','fs'};
    for fi=1:numel(required)
        assert(isfield(R,required{fi}), '%s lacks runResults.%s.',key,required{fi});
    end

    [h,t,l,~,segments,fs] = loadRunSignals(R,key,dataset,opts);
    if isempty(segments)
        warning('pitchControlModel:noFrames','%s contributes no valid segment.',key);
        continue
    end
    windowFrames = max(1,round(opts.changeWindowSec*fs));
    smoothFrames = smoothingWindow(opts.smoothSec,fs);

    for di=1:numel(requestedDelays)
        delayFrames = max(0,round(requestedDelays(di)*fs));
        M = runModels(h,t,l,segments,smoothFrames,windowFrames,delayFrames, ...
            opts.minPairsPerRun);
        runRows{end+1} = modelRow(key,dataset,animal,fs,requestedDelays(di), ... %#ok<AGROW>
            delayFrames/fs,windowFrames/fs,smoothFrames/fs,M);
    end

    M0 = runModels(h,t,l,segments,smoothFrames,windowFrames, ...
        round(opts.primaryCorrectionDelaySec*fs),opts.minPairsPerRun);
    fprintf('  %-52s bLambda=%+5.2f bTorso=%+5.2f | gD=%5.2f | kAlpha=%+6.3f\n', ...
        key,M0.disturbance.bLambda,M0.disturbance.bTorso, ...
        M0.counter.gDisturbance,M0.correction.kAdjusted);
end

byRun = vertcatSafe(runRows);
assert(~isempty(byRun), 'No requested sequence contributed model data.');
changeRmsByRun = primaryChangeRms(byRun,opts.primaryCorrectionDelaySec);
changeRmsSummary = summariseChangeRms(changeRmsByRun);
byAnimal = aggregateByAnimal(byRun);
lagSummary = aggregateLagSummary(byAnimal);
primary = primarySummary(byAnimal,opts.primaryCorrectionDelaySec);

writetable(byRun,fullfile(opts.outputDir,'pitch_control_model_by_run.csv'));
writetable(changeRmsByRun,fullfile(opts.outputDir,'pitch_change_rms_by_run.csv'));
writetable(changeRmsSummary,fullfile(opts.outputDir,'pitch_change_rms_summary.csv'));
writetable(byAnimal,fullfile(opts.outputDir,'pitch_control_model_by_animal.csv'));
writetable(lagSummary,fullfile(opts.outputDir,'pitch_error_correction_lag_summary.csv'));
writetable(primary,fullfile(opts.outputDir,'pitch_control_model_summary.csv'));
plotRepresentativeTimeSeries(opts);
plotContinuousDisturbanceRejection(datasetList,byRun,opts,'run',false);
writeReferenceAdjustedDisturbanceRejection(datasetList,opts);
plotRunLevelGainDistribution(byRun,opts);

A = struct('opts',opts,'byRun',byRun,'changeRmsByRun',changeRmsByRun, ...
    'changeRmsSummary',changeRmsSummary,'byAnimal',byAnimal, ...
    'lagSummary',lagSummary,'primarySummary',primary,'outputDir',opts.outputDir);
save(fullfile(opts.outputDir,'pitch_control_model_analysis.mat'),'A');
writeNotes(A);
end

function [h,t,l,frames,segments,fs] = loadRunSignals(R,key,dataset,opts)
h=R.theta_h_pitch(:); t=R.theta_t_pitch(:); l=R.lambda_h_pitch(:);
mask=logical(R.pitchAnalysisMask(:));
n=min([numel(h),numel(t),numel(l),numel(mask)]);
h=h(1:n);t=t(1:n);l=l(1:n);mask=mask(1:n);
frames=absoluteFrames(R,n);
if opts.excludeInterpolated && isfield(R,'isInterp')
    z=false(n,1);v=logical(R.isInterp(:));z(1:min(n,numel(v)))=v(1:min(n,numel(v)));
    mask=mask&~z;
end
ranges=lookupRanges(key,dataset,opts.frameRanges);
if ~isempty(ranges)
    if strcmp(opts.frameRangeMode,'absFrame'),rv=frames;else,rv=(1:n).';end
    mask=mask&inRanges(rv,ranges);
end
mask=mask&isfinite(h)&isfinite(t)&isfinite(l)&isfinite(frames);
segments=consecutiveSegments(find(mask),frames);
fs=double(R.fs(1));
end

function M = runModels(h,t,l,segments,smoothFrames,windowFrames,delayFrames,minPairs)
nowH=[];nowT=[];nowL=[];
futureAlpha=[];futureQ=[];futureT=[];futureL=[];
for si=1:numel(segments)
    idx=segments{si};
    hs=smoothAngle(h(idx),smoothFrames);
    ts=smoothAngle(t(idx),smoothFrames);
    ls=smoothAngle(l(idx),smoothFrames);
    q=hs-ts;

    if numel(idx)>windowFrames
        a=(1:numel(idx)-windowFrames).';b=a+windowFrames;
        nowH=[nowH;rad2deg(hs(b)-hs(a))]; %#ok<AGROW>
        nowT=[nowT;rad2deg(ts(b)-ts(a))]; %#ok<AGROW>
        nowL=[nowL;rad2deg(ls(b)-ls(a))]; %#ok<AGROW>
    end

    last=numel(idx)-delayFrames-windowFrames;
    if last>=1
        i=(1:last).';a=i+delayFrames;b=a+windowFrames;
        futureAlpha=[futureAlpha;wrapDeg(rad2deg(hs(i)-ls(i)))]; %#ok<AGROW>
        futureQ=[futureQ;rad2deg(q(b)-q(a))]; %#ok<AGROW>
        futureT=[futureT;rad2deg(ts(b)-ts(a))]; %#ok<AGROW>
        futureL=[futureL;rad2deg(ls(b)-ls(a))]; %#ok<AGROW>
    end
end

M=struct();
F=linearFit(nowH,[nowL nowT],minPairs);
M.disturbance=struct('n',F.n,'bLambda',coef(F,2),'bTorso',coef(F,3), ...
    'r2',F.r2,'rmseDeg',F.rmse,'conditionNumber',F.conditionNumber, ...
    'headChangeRmsDeg',rms0(nowH),'torsoChangeRmsDeg',rms0(nowT), ...
    'referenceChangeRmsDeg',rms0(nowL), ...
    'headToTorsoRmsRatio',rms0(nowH)/max(rms0(nowT),eps));
F=linearFit(rad2degDifference(nowH,nowT),[nowT nowL],minPairs);
M.counter=struct('n',F.n,'gDisturbance',-coef(F,2),'gReference',coef(F,3), ...
    'r2',F.r2,'rmseDeg',F.rmse,'conditionNumber',F.conditionNumber);
Fs=linearFit(futureQ,futureAlpha,minPairs);
Fa=linearFit(futureQ,[futureAlpha futureT futureL],minPairs);
M.correction=struct('n',Fa.n,'kSimple',-coef(Fs,2),'simpleR2',Fs.r2, ...
    'kAdjusted',-coef(Fa,2),'gFutureTorso',-coef(Fa,3), ...
    'gFutureReference',coef(Fa,4),'adjustedR2',Fa.r2, ...
    'adjustedRmseDeg',Fa.rmse,'conditionNumber',Fa.conditionNumber, ...
    'opposingFraction',opposingFraction(futureAlpha,futureQ));
end

function dQ = rad2degDifference(dH,dT)
% q=thetaH-thetaT, so the finite change is exactly dH-dT.
dQ=dH-dT;
end

function F = linearFit(y,predictors,minPairs)
y=y(:);X=[ones(numel(y),1),predictors];
good=isfinite(y)&all(isfinite(X),2);y=y(good);X=X(good,:);
F=struct('n',numel(y),'b',nan(size(X,2),1),'r2',NaN,'rmse',NaN, ...
    'conditionNumber',NaN);
if numel(y)<minPairs || size(X,1)<=size(X,2),return;end
F.b=X\y;yhat=X*F.b;res=y-yhat;
F.rmse=sqrt(mean(res.^2));sst=sum((y-mean(y)).^2);
if sst>eps,F.r2=1-sum(res.^2)/sst;end
Z=X(:,2:end);
if isempty(Z)
    F.conditionNumber=1;
else
    Z=(Z-mean(Z,1))./max(std(Z,0,1),eps);
    F.conditionNumber=cond([ones(size(Z,1),1),Z]);
end
end

function value = coef(F,index)
if numel(F.b)>=index,value=F.b(index);else,value=NaN;end
end

function f = opposingFraction(alpha,dq)
good=isfinite(alpha)&isfinite(dq);
if ~any(good),f=NaN;else,f=mean(alpha(good).*dq(good)<0);end
end

function T = modelRow(key,dataset,animal,fs,requestedDelay,actualDelay,window,smooth,M)
T=table(string(key),string(dataset),string(animal),fs,requestedDelay,actualDelay, ...
    window,smooth,M.disturbance.n,M.disturbance.bLambda,M.disturbance.bTorso, ...
    M.disturbance.r2,M.disturbance.rmseDeg,M.disturbance.conditionNumber, ...
    M.disturbance.headChangeRmsDeg,M.disturbance.torsoChangeRmsDeg, ...
    M.disturbance.referenceChangeRmsDeg,M.disturbance.headToTorsoRmsRatio, ...
    M.counter.n,M.counter.gDisturbance,M.counter.gReference,M.counter.r2, ...
    M.counter.rmseDeg,M.counter.conditionNumber,M.correction.n, ...
    M.correction.kSimple,M.correction.simpleR2,M.correction.kAdjusted, ...
    M.correction.gFutureTorso,M.correction.gFutureReference, ...
    M.correction.adjustedR2,M.correction.adjustedRmseDeg, ...
    M.correction.conditionNumber,M.correction.opposingFraction, ...
    'VariableNames',{'datasetKey','dataset','animal','fs','requestedCorrectionDelaySec', ...
    'actualCorrectionDelaySec','changeWindowSec','smoothWindowSec', ...
    'nDisturbancePairs','referenceTrackingGain','torsoTransmissionGain', ...
    'disturbanceModelR2','disturbanceModelRmseDeg','disturbanceConditionNumber', ...
    'headChangeRmsDeg','torsoChangeRmsDeg','referenceChangeRmsDeg', ...
    'headToTorsoRmsRatio', ...
    'nCounterRotationPairs','counterRotationDisturbanceGain','counterRotationReferenceGain', ...
    'counterRotationModelR2','counterRotationModelRmseDeg','counterRotationConditionNumber', ...
    'nCorrectionPairs','errorCorrectionGainSimple','errorCorrectionSimpleR2', ...
    'errorCorrectionGainAdjusted','futureTorsoCompensationGain', ...
    'futureReferenceTrackingGain','errorCorrectionAdjustedR2', ...
    'errorCorrectionAdjustedRmseDeg','errorCorrectionConditionNumber', ...
    'errorOpposingMovementFraction'});
end

function T = primaryChangeRms(byRun,delay)
m=abs(byRun.requestedCorrectionDelaySec-delay)<1e-12;
names={'datasetKey','dataset','animal','nDisturbancePairs','changeWindowSec', ...
    'smoothWindowSec','headChangeRmsDeg','torsoChangeRmsDeg', ...
    'referenceChangeRmsDeg','headToTorsoRmsRatio'};
T=byRun(m,names);
end

function S = summariseChangeRms(T)
S=table(height(T),sum(T.nDisturbancePairs), ...
    mean(T.headChangeRmsDeg,'omitnan'), ...
    mean(T.torsoChangeRmsDeg,'omitnan'), ...
    mean(T.referenceChangeRmsDeg,'omitnan'), ...
    100*mean(T.headToTorsoRmsRatio,'omitnan'), ...
    'VariableNames',{'nRuns','nChangeIntervals','headChangeRmsDeg', ...
    'torsoChangeRmsDeg','referenceChangeRmsDeg','headToTorsoRmsPercent'});
end

function A = aggregateByAnimal(T)
animals=unique(T.animal,'stable');delays=unique(T.requestedCorrectionDelaySec,'stable');rows={};
metrics=metricNames();
for di=1:numel(delays)
    for ai=1:numel(animals)
        m=T.animal==animals(ai)&T.requestedCorrectionDelaySec==delays(di);
        if ~any(m),continue;end
        values=cell(1,numel(metrics));
        for j=1:numel(metrics),values{j}=mean(T.(metrics{j})(m),'omitnan');end
        rows{end+1}=cell2table([{animals(ai),delays(di),mean(T.actualCorrectionDelaySec(m)), ... %#ok<AGROW>
            sum(m)},values],'VariableNames',[{'animal','requestedCorrectionDelaySec', ...
            'meanActualCorrectionDelaySec','nRuns'},metrics]);
    end
end
A=vertcatSafe(rows);
end

function names = metricNames()
names={'referenceTrackingGain','torsoTransmissionGain','disturbanceModelR2', ...
    'counterRotationDisturbanceGain','counterRotationReferenceGain', ...
    'counterRotationModelR2','errorCorrectionGainSimple', ...
    'errorCorrectionGainAdjusted','futureTorsoCompensationGain', ...
    'futureReferenceTrackingGain','errorCorrectionAdjustedR2', ...
    'errorOpposingMovementFraction'};
end

function S = aggregateLagSummary(A)
delays=unique(A.requestedCorrectionDelaySec,'stable');rows={};
for i=1:numel(delays)
    m=A.requestedCorrectionDelaySec==delays(i);
    ks=A.errorCorrectionGainSimple(m);ka=A.errorCorrectionGainAdjusted(m);
    rows{end+1}=table(delays(i),mean(A.meanActualCorrectionDelaySec(m)),sum(m), ... %#ok<AGROW>
        mean(ks,'omitnan'),se(ks),mean(ka,'omitnan'),se(ka), ...
        mean(A.errorOpposingMovementFraction(m),'omitnan'), ...
        'VariableNames',{'requestedCorrectionDelaySec','meanActualCorrectionDelaySec', ...
        'nAnimals','errorCorrectionGainSimple','errorCorrectionGainSimpleSE', ...
        'errorCorrectionGainAdjusted','errorCorrectionGainAdjustedSE', ...
        'errorOpposingMovementFraction'});
end
S=vertcatSafe(rows);
end

function S = primarySummary(A,delay)
m=abs(A.requestedCorrectionDelaySec-delay)<1e-12;
P=A(m,:);names=metricNames();metric=strings(0,1);meanValue=[];seValue=[];n=[];
for i=1:numel(names)
    v=P.(names{i});metric(end+1,1)=string(names{i}); %#ok<AGROW>
    meanValue(end+1,1)=mean(v,'omitnan');seValue(end+1,1)=se(v); %#ok<AGROW>
    n(end+1,1)=nnz(isfinite(v)); %#ok<AGROW>
end
S=table(repmat(delay,numel(names),1),metric,n,meanValue,seValue, ...
    'VariableNames',{'correctionDelaySec','metric','nAnimals','animalEqualMean','animalSE'});
end

function plotRepresentativeTimeSeries(opts)
% Panel B: one pre-designated long, continuous sequence for visual context.
key=normaliseKey(opts.representativeDataset);
[group,dataset]=splitDatasetKey(key);
path=fullfile(opts.dataRoot,group,dataset,opts.resultFileName);
assert(isfile(path),'Missing representative sequence result: %s',path);
S=load(path,'runResults');R=S.runResults;
[h,t,l,frames,segments,fs]=loadRunSignals(R,key,dataset,opts);
assert(~isempty(segments),'Representative sequence has no valid pitch segment: %s',key);
[~,which]=max(cellfun(@numel,segments));idx=segments{which};
w=smoothingWindow(opts.representativeSmoothSec,fs);
h=rad2deg(smoothAngle(h(idx),w));t=rad2deg(smoothAngle(t(idx),w));
l=rad2deg(smoothAngle(l(idx),w));q=h-t;alpha=wrapDeg(h-l);
timeSec=(frames(idx)-frames(idx(1)))/fs;

T=table(repmat(string(key),numel(idx),1),frames(idx),timeSec,h,t,l,q,alpha, ...
    'VariableNames',{'datasetKey','absFrame','timeSec','thetaHdeg','thetaTdeg', ...
    'lambdaDeg','qDeg','alphaDeg'});
writetable(T,fullfile(opts.outputDir,'pitch_representative_control_timeseries.csv'));

visible=onoff(opts.showFigures);f=figure('Color','w','Visible',visible, ...
    'Position',[100 100 1120 500]);ax=axes(f);hold(ax,'on');
colors=struct('torso',[.88 .30 .07],'relative',[.20 .62 .35], ...
    'head',[.05 .43 .72],'reference',[.12 .12 .12],'error',[.55 .25 .65]);
pT=plot(ax,timeSec,t,'-','Color',colors.torso,'LineWidth',2.0);
pQ=plot(ax,timeSec,q,'-','Color',colors.relative,'LineWidth',2.0);
pH=plot(ax,timeSec,h,'-','Color',colors.head,'LineWidth',2.25);
pL=plot(ax,timeSec,l,'--','Color',colors.reference,'LineWidth',1.8);
pA=plot(ax,timeSec,alpha,':','Color',colors.error,'LineWidth',2.0);
yline(ax,0,':','Color',[.55 .55 .55],'HandleVisibility','off');
grid(ax,'on');box(ax,'off');xlabel(ax,'Time from start of valid segment (s)');
ylabel(ax,'Pitch angle (deg)');
title(ax,sprintf('Representative pursuit sequence: %s',dataset), ...
    'Interpreter','none');
legend(ax,[pT pQ pH pL pA],{'Torso pitch, \theta_T','Head-on-torso pitch, q', ...
    'World-frame head pitch, \theta_H','Head-lure-line elevation, \lambda','Pitch error, \alpha'}, ...
    'Interpreter','tex','Location','northeast');
allValues=[t;q;h;l;alpha];lim=ceil(max(abs(allValues))/5)*5;ylim(ax,[-lim lim]);
saveFigurePair(f,fullfile(opts.outputDir,'pitch_representative_control_timeseries'),opts);
end

function plotContinuousDisturbanceRejection(datasetList,A,opts,weighting,showTitle)
% Panel C sign convention: y=Delta q and x=-Delta theta_T. A unit positive
% slope therefore means complete counter-rotation, exactly as in the plan.
rows={};runs=struct('key',{},'animal',{},'x',{},'y',{});
for ri=1:numel(datasetList)
    key=normaliseKey(datasetList{ri});[group,dataset]=splitDatasetKey(key);
    animal=inferAnimal(dataset);path=fullfile(opts.dataRoot,group,dataset,opts.resultFileName);
    S=load(path,'runResults');R=S.runResults;
    [h,t,~,~,segments,fs]=loadRunSignals(R,key,dataset,opts);
    w=smoothingWindow(opts.smoothSec,fs);lag=max(1,round(opts.changeWindowSec*fs));
    runX=[];runY=[];
    for si=1:numel(segments)
        idx=segments{si};hs=smoothAngle(h(idx),w);ts=smoothAngle(t(idx),w);q=hs-ts;
        if numel(idx)<=lag,continue;end
        a=(1:numel(idx)-lag).';b=a+lag;
        runX=[runX;-rad2deg(ts(b)-ts(a))]; %#ok<AGROW>
        runY=[runY; rad2deg(q(b)-q(a))]; %#ok<AGROW>
    end
    good=isfinite(runX)&isfinite(runY);runX=runX(good);runY=runY(good);
    rows{end+1}=table(repmat(string(key),numel(runX),1), ... %#ok<AGROW>
        repmat(string(animal),numel(runX),1),runX,runY, ...
        'VariableNames',{'datasetKey','animal','opposingTorsoChangeDeg','qChangeDeg'});
    runs(end+1)=struct('key',string(key),'animal',string(animal), ... %#ok<AGROW>
        'x',runX,'y',runY);
end
writetable(vertcatSafe(rows), ...
    fullfile(opts.outputDir,'pitch_continuous_disturbance_rejection.csv'));

m=abs(A.requestedCorrectionDelaySec-opts.primaryCorrectionDelaySec)<1e-12;
gain=A.counterRotationDisturbanceGain(m);mu=mean(gain,'omitnan');s=se(gain);
allX=vertcat(runs.x);allY=vertcat(runs.y);
lim=max(10,ceil(percentileLocal(abs([allX;allY]),99)/5)*5);edges=linspace(-lim,lim,51);
centres=(edges(1:end-1)+edges(2:end))/2;
N=weightedDensity(runs,edges,weighting);
green=[.20 .62 .35];visible=onoff(opts.showFigures);
f=figure('Color','w','Visible',visible,'Position',[100 100 700 650]);ax=axes(f);hold(ax,'on');
im=imagesc(ax,centres,centres,N.');set(im,'AlphaData',.10+.90*(N.'>0));
set(ax,'YDir','normal');tmap=linspace(0,1,256).';colormap(ax,(1-tmap).*ones(256,3)+tmap.*green);
pIdentity=plot(ax,[-lim lim],[-lim lim],'--','Color',[.15 .15 .15],'LineWidth',1.5);
pGain=plot(ax,[-lim lim],mu*[-lim lim],'-','Color',green,'LineWidth',2.4);
xline(ax,0,':','Color',[.55 .55 .55],'HandleVisibility','off');
yline(ax,0,':','Color',[.55 .55 .55],'HandleVisibility','off');
axis(ax,'square');xlim(ax,[-lim lim]);ylim(ax,[-lim lim]);grid(ax,'on');box(ax,'off');
xlabel(ax,'Opposing torso-pitch change, -\Delta\theta_T (deg)');
ylabel(ax,'Head-on-torso pitch change, \Deltaq (deg)');
if strcmp(weighting,'run')
    weightingLabel='Run-equal';base='pitch_continuous_disturbance_rejection_run_equal';
else
    weightingLabel='Animal-equal';base='pitch_continuous_disturbance_rejection';
end
if showTitle
    title(ax,sprintf('%s continuous torso-disturbance rejection',weightingLabel));
end
legend(ax,[pIdentity pGain],{'Complete cancellation (g_D = 1)', ...
    sprintf('%s g_D = %.2f \\pm %.2f SE',weightingLabel,mu,s)}, ...
    'Location','northwest','Box','on','Color','white', ...
    'EdgeColor',[.15 .15 .15],'Interpreter','tex');
saveFigurePair(f,fullfile(opts.outputDir,base),opts);
end

function N=weightedDensity(runs,edges,weighting)
% Each run first contributes a unit-area density. Units are then averaged:
% runs for run-equal, or run-equal within animal followed by animal-equal.
nBins=numel(edges)-1;N=zeros(nBins,nBins);
if strcmp(weighting,'run')
    units=string({runs.key});
else
    units=string({runs.animal});
end
units=unique(units,'stable');nUsed=0;
for ui=1:numel(units)
    if strcmp(weighting,'run')
        idx=find(string({runs.key})==units(ui));
    else
        idx=find(string({runs.animal})==units(ui));
    end
    unitN=zeros(nBins,nBins);nRuns=0;
    for ri=idx
        H=histcounts2(runs(ri).x,runs(ri).y,edges,edges);
        if sum(H,'all')==0,continue;end
        unitN=unitN+H/sum(H,'all');nRuns=nRuns+1;
    end
    if nRuns==0,continue;end
    N=N+unitN/nRuns;nUsed=nUsed+1;
end
if nUsed>0,N=N/nUsed;end
end

function writeReferenceAdjustedDisturbanceRejection(datasetList,opts)
% Retain the diagnostic table without exporting a non-manuscript figure.
% Within each run, both y=Delta q and x=-Delta theta_T are residualised
% against concurrent Delta lambda. By the Frisch--Waugh--Lovell theorem,
% the resulting slope equals the adjusted g_D from the full model.
supplementaryDir=fullfile(opts.outputDir,'supplementary');
if ~exist(supplementaryDir,'dir'),mkdir(supplementaryDir);end
rows={};runs=struct('key',{},'animal',{},'x',{},'y',{});
for ri=1:numel(datasetList)
    key=normaliseKey(datasetList{ri});[group,dataset]=splitDatasetKey(key);
    animal=inferAnimal(dataset);path=fullfile(opts.dataRoot,group,dataset,opts.resultFileName);
    S=load(path,'runResults');R=S.runResults;
    [h,t,l,~,segments,fs]=loadRunSignals(R,key,dataset,opts);
    w=smoothingWindow(opts.smoothSec,fs);lag=max(1,round(opts.changeWindowSec*fs));
    runX=[];runY=[];runL=[];
    for si=1:numel(segments)
        idx=segments{si};hs=smoothAngle(h(idx),w);ts=smoothAngle(t(idx),w);
        ls=smoothAngle(l(idx),w);q=hs-ts;
        if numel(idx)<=lag,continue;end
        a=(1:numel(idx)-lag).';b=a+lag;
        runX=[runX;-rad2deg(ts(b)-ts(a))]; %#ok<AGROW>
        runY=[runY; rad2deg(q(b)-q(a))]; %#ok<AGROW>
        runL=[runL; rad2deg(ls(b)-ls(a))]; %#ok<AGROW>
    end
    good=isfinite(runX)&isfinite(runY)&isfinite(runL);
    runX=runX(good);runY=runY(good);runL=runL(good);
    Z=[ones(numel(runL),1) runL];
    residualX=runX-Z*(Z\runX);
    residualY=runY-Z*(Z\runY);
    rows{end+1}=table(repmat(string(key),numel(runX),1), ... %#ok<AGROW>
        repmat(string(animal),numel(runX),1),runX,runY,runL,residualX,residualY, ...
        'VariableNames',{'datasetKey','animal','opposingTorsoChangeDeg', ...
        'qChangeDeg','referenceChangeDeg','opposingTorsoResidualDeg','qResidualDeg'});
    runs(end+1)=struct('key',string(key),'animal',string(animal), ... %#ok<AGROW>
        'x',residualX,'y',residualY);
end
writetable(vertcatSafe(rows),fullfile(supplementaryDir, ...
    'pitch_reference_adjusted_disturbance_rejection.csv'));

end

function plotRunLevelGainDistribution(A,opts)
% Compact run-level torso-to-head transmission gain, nested within animal.
m=abs(A.requestedCorrectionDelaySec-opts.primaryCorrectionDelaySec)<1e-12;
P=A(m,:);assert(~isempty(P),'Primary correction delay is absent from run table.');
animals=unique(P.animal,'stable');nAnimals=numel(animals);
g=P.torsoTransmissionGain;
cols=lines(nAnimals);visible=onoff(opts.showFigures);
f=figure('Color','w','Visible',visible,'Position',[100 100 800 560]);
ax=axes(f,'Position',[0.13 0.18 0.75 0.72]);hold(ax,'on');
set(ax,'FontSize',12);
for ai=1:nAnimals
    idx=strcmp(P.animal,animals{ai});v=g(idx);n=sum(idx);
    offsets=linspace(-.17,.17,n);
    x=ai+offsets;
    animalMean=mean(v,'omitnan');
    for ri=1:n
        h=plot(ax,[x(ri) ai],[v(ri) animalMean],'-', ...
            'Color',.32*cols(ai,:)+.68*[1 1 1], ...
            'LineWidth',.8);
        if ai==1 && ri==1,hLink=h;end
    end
    h=plot(ax,x,v,'o','MarkerSize',5.5,'MarkerFaceColor',cols(ai,:), ...
        'MarkerEdgeColor','w','LineWidth',.7);
    if ai==1,hRun=h;end
    h=plot(ax,ai,animalMean,'o','MarkerSize',11,'MarkerFaceColor',cols(ai,:), ...
        'MarkerEdgeColor',[.1 .1 .1],'LineWidth',1.2);
    if ai==1,hAnimal=h;end
end
animalMeans=nan(nAnimals,1);
for ai=1:nAnimals
    animalMeans(ai)=mean(g(strcmp(P.animal,animals{ai})),'omitnan');
end
populationMean=mean(animalMeans,'omitnan');populationSE=se(animalMeans);
populationX=nAnimals+1.15;
hPopulation=errorbar(ax,populationX,populationMean,populationSE,'d','MarkerSize',10, ...
    'MarkerFaceColor','w','MarkerEdgeColor',[.1 .1 .1],'Color',[.1 .1 .1], ...
    'LineWidth',1.8,'CapSize',10);
yline(ax,0,'--','g_T = 0: complete rejection','Color',[.15 .15 .15],'LineWidth',1.2, ...
    'LabelHorizontalAlignment','left');
yline(ax,1,'--','g_T = 1: complete transmission','Color',[.45 .45 .45], ...
    'LineWidth',1.2,'LabelHorizontalAlignment','left');
xlim(ax,[.45 populationX+.95]);ylim(ax,coefficientLimits([g;animalMeans;0;1], ...
    0,populationMean,populationSE));
cheetahLabels=arrayfun(@(i)sprintf('Cheetah %d',i),1:nAnimals, ...
    'UniformOutput',false);
set(ax,'XTick',[1:nAnimals populationX], ...
    'XTickLabel',[cheetahLabels {'Pop.'}], ...
    'TickLabelInterpreter','none','YGrid','on');
xtickangle(ax,0);box(ax,'off');
ylabel(ax,'g_T','Interpreter','tex');
h=legend(ax,[hLink hRun hAnimal hPopulation], ...
    {'Run-to-animal mean','Individual run','Animal mean', ...
    sprintf('Population: g_T = %.2f \\pm %.2f SE',populationMean,populationSE)}, ...
    'Location','northeast','Box','on','Color','white','Interpreter','tex');
h.FontSize=9;
h.ItemTokenSize=[18 8];
h.Units='normalized';
drawnow;
legendPosition=h.Position;
legendPosition(2)=0.57;
h.Location='none';
h.Position=legendPosition;
h.EdgeColor=[.15 .15 .15];
h.AutoUpdate='off';
saveFigurePair(f,fullfile(opts.outputDir,'pitch_population_control_gains_run_level'),opts);
end

function lim=coefficientLimits(v,reference,mu,s)
x=[v(:);reference;mu-s;mu+s];x=x(isfinite(x));
if isempty(x),lim=[-1 1];return;end
span=max(x)-min(x);pad=max(.15,.12*max(span,eps));lim=[min(x)-pad max(x)+pad];
end

function writeNotes(A)
path=fullfile(A.outputDir,'PITCH_MECHANISM_RESULTS.md');fid=fopen(path,'w');
assert(fid>0,'Cannot create %s.',path);c=onCleanup(@()fclose(fid));
S=A.primarySummary;
fprintf(fid,'# Pitch mechanism analysis\n\n');
fprintf(fid,'All coefficients are first fitted per pursuit sequence, averaged run-equal within animal, and then summarised across animals. Changes use a %.0f-ms window after %.0f-ms smoothing.\n\n', ...
    1000*A.opts.changeWindowSec,1000*A.opts.smoothSec);
fprintf(fid,'## Sign convention\n\n');
fprintf(fid,'- `q = theta_H - theta_T` (head-on-torso pitch).\n');
fprintf(fid,'- `alpha = theta_H - lambda` (measured pitch error); controller error is `e = -alpha`.\n');
fprintf(fid,'- Counter-rotation is plotted as `Delta q` versus `-Delta theta_T`, so complete cancellation has slope `g_D = +1`.\n');
fprintf(fid,'- Error correction is `Delta q(t+tau) = -k_alpha alpha(t)`, so corrective responses have `k_alpha > 0`.\n\n');
fprintf(fid,'## Figure panels\n\n');
fprintf(fid,'- `pitch_representative_control_timeseries.png` (theta_T, q, theta_H, head-lure-line elevation lambda, and alpha in the full named representative sequence).\n');
fprintf(fid,'- `pitch_continuous_disturbance_rejection.png` (animal-equal Delta q versus -Delta theta_T); `pitch_continuous_disturbance_rejection_run_equal.png` is the run-equal companion.\n');
fprintf(fid,'- `pitch_population_control_gains_run_level.png` (per-run torso-to-head transmission gains).\n\n');
fprintf(fid,'- Supplementary: `supplementary/pitch_reference_adjusted_disturbance_rejection.png` visualises the partial g_D relationship after residualising both axes within run against Delta lambda.\n\n');
fprintf(fid,'## Primary coefficient results\n\n');
writeMetric(fid,S,'torsoTransmissionGain','Torso transmission gain, b_T',0);
writeMetric(fid,S,'counterRotationDisturbanceGain','Counter-rotation gain, g_D',1);
writeMetric(fid,S,'referenceTrackingGain','Reference tracking gain, b_lambda',1);
writeMetric(fid,S,'errorCorrectionGainAdjusted','Adjusted delayed error-correction gain, k_alpha',NaN);
R=A.changeRmsSummary;
fprintf(fid,['\nRun-equal %.0f-ms change RMS: head %.3f deg, torso %.3f deg, ' ...
    'reference %.3f deg; mean within-run head/torso RMS ratio %.3f%% ' ...
    '(n=%d runs, %d intervals).\n'],1000*A.opts.changeWindowSec, ...
    R.headChangeRmsDeg,R.torsoChangeRmsDeg,R.referenceChangeRmsDeg, ...
    R.headToTorsoRmsPercent,R.nRuns,R.nChangeIntervals);
fprintf(fid,['\nBecause q=theta_H-theta_T by definition in this planar analysis, the two ' ...
    'finite-change models obey g_D=1-b_T and g_lambda=b_lambda. The ' ...
    'counter-rotation coefficients describe the kinematic implementation of ' ...
    'disturbance rejection; they are not statistically independent evidence.\n']);
fprintf(fid,'\nThe error-correction response window begins %.0f ms after the error measurement. Its adjusted model controls torso and reference changes during that future window. Positive k_alpha has the predicted corrective sign. This observational lagged association does not by itself identify neural feedback or exclude predictive/mechanical control.\n', ...
    1000*A.opts.primaryCorrectionDelaySec);
end

function writeMetric(fid,S,name,label,prediction)
m=S.metric==name;
if any(m)
    if isfinite(prediction),suffix=sprintf(' (prediction %.1f)',prediction);else,suffix='';end
    fprintf(fid,'- %s: %.3f +/- %.3f animal SE, n=%d%s.\n',label, ...
        S.animalEqualMean(m),S.animalSE(m),S.nAnimals(m),suffix);
end
end

function value = se(x)
x=x(isfinite(x));if numel(x)<2,value=NaN;else,value=std(x)/sqrt(numel(x));end
end
function y=rms0(x)
x=x(isfinite(x));if isempty(x),y=NaN;else,y=sqrt(mean(x.^2));end
end
function q=percentileLocal(x,p)
x=sort(x(isfinite(x)));if isempty(x),q=NaN;return;end
if numel(x)==1,q=x;return;end
k=1+(numel(x)-1)*p/100;lo=floor(k);hi=ceil(k);q=x(lo)+(k-lo)*(x(hi)-x(lo));
end
function y=wrapDeg(x),y=mod(x+180,360)-180;end
function x=smoothAngle(x,w)
x=unwrap(x(:));
if w>1&&numel(x)>=3,w=min(w,numel(x));if mod(w,2)==0,w=w-1;end
    if w>=3,x=movmedian(x,w,'omitmissing','Endpoints','shrink');end
end
end
function w=smoothingWindow(sec,fs)
if sec<=0,w=1;return;end
w=max(3,round(sec*fs));if mod(w,2)==0,w=w+1;end
end
function segments=consecutiveSegments(indices,frames)
segments={};if isempty(indices),return;end
cut=[true;diff(indices(:))~=1|diff(frames(indices))~=1];s=find(cut);e=[s(2:end)-1;numel(indices)];
for i=1:numel(s),segments{end+1}=indices(s(i):e(i));end %#ok<AGROW>
end
function frames=absoluteFrames(R,n)
if isfield(R,'commonFrames')&&numel(R.commonFrames)>=n
    frames=double(R.commonFrames(1:n));frames=frames(:);
elseif isfield(R,'absFrames')&&numel(R.absFrames)>=n
    frames=double(R.absFrames(1:n));frames=frames(:);
else
    frames=(1:n).';
end
end
function ranges=lookupRanges(key,dataset,entries)
ranges=[];if isempty(entries),return;end
for i=1:size(entries,1),candidate=normaliseKey(entries{i,1});
    if strcmp(candidate,key)||strcmp(candidate,dataset),r=double(entries{i,2});if isvector(r),r=reshape(r,1,2);end;ranges=[ranges;sort(r,2)];end %#ok<AGROW>
end
end
function tf=inRanges(v,ranges)
tf=false(size(v));for i=1:size(ranges,1),tf=tf|(v>=ranges(i,1)&v<=ranges(i,2));end
end
function animal=inferAnimal(dataset)
token=regexp(char(dataset),'^(.*?)(?:Top|Bottom|Run|Flick)','tokens','once');
if isempty(token),animal=string(dataset);else,animal=string(token{1});end
end
function list=normaliseList(list)
if isstring(list),list=cellstr(list(:));elseif ischar(list),list={list};end;list=list(:);
end
function key=normaliseKey(key)
key=strrep(char(key),'\','/');key=regexprep(key,'^/+|/+$','');
end
function [group,dataset]=splitDatasetKey(key)
p=strsplit(key,'/');assert(numel(p)>=2,'Expected group/dataset: %s',key);group=p{end-1};dataset=p{end};
end
function T=vertcatSafe(rows)
if isempty(rows),T=table();else,T=vertcat(rows{:});end
end
function opts=applyDefaults(opts,defaults)
f=fieldnames(defaults);for i=1:numel(f),if ~isfield(opts,f{i})||isempty(opts.(f{i})),opts.(f{i})=defaults.(f{i});end;end
end
function s=onoff(tf),if tf,s='on';else,s='off';end,end
function shade(ax,x,lo,hi,c)
good=isfinite(x)&isfinite(lo)&isfinite(hi);if nnz(good)<2,return;end
fill(ax,[x(good);flipud(x(good))],[lo(good);flipud(hi(good))],c, ...
    'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off');
end
function saveFigurePair(f,base,opts)
if opts.saveFigures,exportgraphics(f,[base '.png'],'Resolution',opts.figureResolution); ...
        exportgraphics(f,[base '.pdf'],'ContentType','vector');end
if ~opts.showFigures,close(f);end
end
