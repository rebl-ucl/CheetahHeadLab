function R=headTorsoHandoff(opts)
% Head/torso separation, range closure and exploratory lagged coupling.
% Runs are the independent units. Coupling is descriptive and does not
% establish that head orientation causes subsequent steering.

arguments
    opts struct=struct()
end

thisDir=fileparts(mfilename('fullpath'));
pursuitDir=fileparts(thisDir);repoRoot=fileparts(thisDir);
opts=applyDefaults(opts,struct( ...
    'dataRoot',fullfile(repoRoot,'outputs','analysis_inputs','yaw_tracking'), ...
    'outputDir',fullfile(repoRoot,'outputs','yaw','head_torso'), ...
    'nMonteCarlo',300, ...
    'tauGridSec',[0.050 0.075 0.100 0.150], ...
    'minPairs',20, ...
    'nSurrogate',2000, ...
    'velocityWindowSec',0.100, ...
    'missVelocityWindowSec',0.200));
if ~exist(opts.outputDir,'dir'),mkdir(opts.outputDir);end

names={'julesFlick1_20190309','julesTopRun11_20170829', ...
    'julesTopRun1_20171209','lilyRun_20190309','phantomRun_20190303'};
labels={'Jules Flick 1','Jules Top Run 11','Jules Top Run 1', ...
    'Lily Run','Phantom Run'};
animals={'Jules','Jules','Jules','Lily','Phantom'};
nRun=numel(names);nDraw=opts.nMonteCarlo;nTau=numel(opts.tauGridSec);

fprintf('\n==== headTorsoHandoff ====\n');
fprintf('Runs: %d | joint reconstruction draws per run: %d | tau grid: %s ms\n', ...
    nRun,nDraw,mat2str(1000*opts.tauGridSec));

D=loadInputs(opts,names);

% ---- nominal reconstruction -------------------------------------------------
for j=1:nRun
    tmp=computeRunMetrics(D(j).nominalPose,D(j).nominalLure,D(j).fs,D(j),opts);
    if j==1,nominal=repmat(tmp,nRun,1);else,nominal(j)=tmp;end
end

% ---- Monte Carlo propagation ------------------------------------------------
scalarFlds={'medAbsTravelMinusTorsoDeg','medAbsTravelMinusHeadDeg', ...
    'travelAlignmentGainDeg','rangeStartM','rangeEndM','rangeChangeM', ...
    'medClosingSpeedMps','fracFramesClosing','medMissDistanceM', ...
    'missDistanceChangeM','medTimeToClosestSec','fracFramesApproaching', ...
    'medAbsEtaDeg','medAbsDeltaHeadDeg', ...
    'medAbsDeltaTorsoDeg','deltaHeadMeanDeg','deltaHeadRmsDeg', ...
    'deltaTorsoMeanDeg','deltaTorsoRmsDeg','lockRatioTorsoOverHead', ...
    'nullConstWorldHeadingTorsoDeg','nullConstWorldHeadingHeadDeg', ...
    'nullFollowLosTorsoDeg','nullFollowLosHeadDeg', ...
    'losLockGainTorso','losLockGainHead','effectiveNTorso'};
vecFlds={'pcHeadTurnToTorso','pcTorsoTurnToHead','pcEtaToTorsoTurn', ...
    'pcHeadLeadToTorsoTurn'};
mc=struct();
for f=1:numel(scalarFlds),mc.(scalarFlds{f})=nan(nDraw,nRun);end
for f=1:numel(vecFlds),mc.(vecFlds{f})=nan(nDraw,nRun,nTau);end
for k=1:nDraw
    for j=1:nRun
        idx=D(j).successfulDrawIndex(k);
        pose=double(squeeze(D(j).poseAll(idx,:,:)));
        lure=double(squeeze(D(j).lureAll(idx,:,:)));
        m=computeRunMetrics(pose,lure,D(j).fs,D(j),opts);
        for f=1:numel(scalarFlds),mc.(scalarFlds{f})(k,j)=m.(scalarFlds{f});end
        for f=1:numel(vecFlds),mc.(vecFlds{f})(k,j,:)=m.(vecFlds{f});end
    end
    if mod(k,50)==0||k==nDraw,fprintf('  propagated %d/%d joint draws\n',k,nDraw);end
end

% ---- circular-shift surrogate null (nominal reconstruction) -----------------
% Preserves the autocorrelation of both series and destroys only their temporal
% alignment, so it measures how large |r| gets by chance in a run this short.
surro=makeSurrogateTable(D,nominal,names,labels,opts);

% ---- summary tables --------------------------------------------------------
runSummary=makeRunSummary(nominal,mc,names,labels,animals,scalarFlds);
couplingTau=makeCouplingTable(nominal,mc,names,labels,animals,vecFlds,opts);
betweenRun=makeBetweenRunSigns(mc,vecFlds,opts);
byFrame=makeByFrameTable(nominal,names,animals);

writetable(runSummary,fullfile(opts.outputDir,'run_handoff_summary.csv'));
writetable(couplingTau,fullfile(opts.outputDir,'coupling_by_tau.csv'));
writetable(betweenRun,fullfile(opts.outputDir,'between_run_sign_tests.csv'));
writetable(surro,fullfile(opts.outputDir,'surrogate_null.csv'));
writetable(byFrame,fullfile(opts.outputDir,'nominal_by_frame.csv'));

R=struct('opts',opts,'runNames',{names},'runLabels',{labels}, ...
    'animals',{animals},'nominal',nominal,'mc',mc,'runSummary',runSummary, ...
    'couplingTau',couplingTau,'betweenRun',betweenRun,'surrogate',surro, ...
    'byFrame',byFrame);
save(fullfile(opts.outputDir,'head_torso_handoff.mat'),'R','-v7.3');
fprintf('Saved head/torso handoff results to %s\n',opts.outputDir);
printHeadline(runSummary,couplingTau,betweenRun,opts);
end

% =========================================================================
function D=loadInputs(opts,names)
D=repmat(struct(),numel(names),1);
for j=1:numel(names)
    name=names{j};folder=fullfile(opts.dataRoot,name);
    mcPath=fullfile(folder,'MonteCarloResults',['mcMain_' name '.mat']);
    assert(isfile(mcPath),'Missing Monte Carlo input: %s',mcPath);
    S=load(mcPath,'poseAll','lureXYZAll','nominalPose','nominalLureXYZ', ...
        'sampleSucceeded','absFrames');
    A=load(fullfile(folder,'lureAnalysis_results_trajectory.mat'),'runResults');
    A=A.runResults;frames=double(S.absFrames(:));
    if isfield(A,'yawAnalysisMask')&&numel(A.yawAnalysisMask)==numel(frames)
        mask=logical(A.yawAnalysisMask(:));
    else
        mask=true(numel(frames),1);
    end
    ok=find(logical(S.sampleSucceeded(:)));
    assert(numel(ok)>=opts.nMonteCarlo, ...
        '%s has only %d successful joint draws; %d requested.', ...
        name,numel(ok),opts.nMonteCarlo);
    D(j).name=name;D(j).frames=frames;D(j).fs=double(A.fs(1));
    D(j).analysisMask=mask;
    D(j).nominalPose=double(S.nominalPose);D(j).nominalLure=double(S.nominalLureXYZ);
    D(j).poseAll=S.poseAll;D(j).lureAll=S.lureXYZAll;
    D(j).successfulDrawIndex=ok(1:opts.nMonteCarlo);
end
end

% =========================================================================
function m=computeRunMetrics(pose,lure,fs,Drun,opts)
% All angles in radians internally; every reported quantity in degrees.
N=size(pose,1);t=(Drun.frames-Drun.frames(1))/fs;
eye=pose(:,4:6);torso=pose(:,10:12);

psiHead =forwardYaw(pose(:,1),pose(:,2),pose(:,3));   % head  world yaw
psiTorso=forwardYaw(pose(:,7),pose(:,8),pose(:,9));   % torso world yaw

lamHead =atan2(lure(:,2)-eye(:,2)  ,lure(:,1)-eye(:,1));    % eye  ->lure bearing
lamTorso=atan2(lure(:,2)-torso(:,2),lure(:,1)-torso(:,1));  % torso->lure bearing

eta       =wrapPi(psiHead-psiTorso);   % gaze offset on the body
deltaHead =wrapPi(psiHead-lamHead);    % head  angular lead on its LOS
deltaTorso=wrapPi(psiTorso-lamTorso);  % torso angular offset on its LOS

uHead=unwrap(psiHead);uTorso=unwrap(psiTorso);
uLamHead=unwrap(lamHead);uLamTorso=unwrap(lamTorso);

% travel heading of the torso origin (causal trailing fit; premise check only)
vx=trailingVelocity(torso(:,1),t,fs,opts.velocityWindowSec);
vy=trailingVelocity(torso(:,2),t,fs,opts.velocityWindowSec);
travel=atan2(vy,vx);

valid=Drun.analysisMask(:)&isfinite(psiHead)&isfinite(psiTorso)& ...
    isfinite(lamHead)&isfinite(lamTorso);

m=struct();m.N=nnz(valid);

% ---- premise: torso is the locomotor axis, head is the sensor --------------
% Both segments are compared to the SAME travel heading, so any error in the
% travel estimate is shared and the paired difference is the robust quantity.
gTrav=valid&isfinite(travel);
travelMinusTorso=wrapPi(travel-psiTorso);
travelMinusHead =wrapPi(travel-psiHead);
m.medAbsTravelMinusTorsoDeg=rad2deg(median(abs(travelMinusTorso(gTrav))));
m.medAbsTravelMinusHeadDeg =rad2deg(median(abs(travelMinusHead(gTrav))));
% >0 means the torso is the axis closer to where the animal is actually going.
m.travelAlignmentGainDeg=m.medAbsTravelMinusHeadDeg-m.medAbsTravelMinusTorsoDeg;
m.medAbsEtaDeg       =rad2deg(median(abs(eta(valid))));
m.medAbsDeltaHeadDeg =rad2deg(median(abs(deltaHead(valid))));
m.medAbsDeltaTorsoDeg=rad2deg(median(abs(deltaTorso(valid))));

% ---- segment LOS lock: mean offset and RMS scatter about it ----------------
m.deltaHeadMeanDeg =rad2deg(circMean(deltaHead(valid)));
m.deltaHeadRmsDeg  =rmsAboutCircMean(deltaHead(valid));
m.deltaTorsoMeanDeg=rad2deg(circMean(deltaTorso(valid)));
m.deltaTorsoRmsDeg =rmsAboutCircMean(deltaTorso(valid));
m.lockRatioTorsoOverHead=m.deltaTorsoRmsDeg/max(m.deltaHeadRmsDeg,eps);

% Trivial-fit nulls.  "Constant world heading" is what the residual would be if
% the segment never turned (all variation left in the LOS); "follow LOS" is what
% it would be if the segment tracked the LOS with no offset dynamics of its own.
m.nullConstWorldHeadingTorsoDeg=rmsAboutCircMean(lamTorso(valid));
m.nullConstWorldHeadingHeadDeg =rmsAboutCircMean(lamHead(valid));
m.nullFollowLosTorsoDeg=rmsAboutCircMean(psiTorso(valid));
m.nullFollowLosHeadDeg =rmsAboutCircMean(psiHead(valid));
% Gain of the LOS-referenced description over simply standing still in world
% yaw.  >1 means "constant offset on the LOS" genuinely beats "constant world
% heading"; ~1 or below means the constant offset is absorbing a geometry that
% barely moved, and the lock claim is not earned.
m.losLockGainTorso=m.nullConstWorldHeadingTorsoDeg/max(m.deltaTorsoRmsDeg,eps);
m.losLockGainHead =m.nullConstWorldHeadingHeadDeg /max(m.deltaHeadRmsDeg ,eps);

% ---- is the locomotion functionally aimed at the lure? ---------------------
% Three progressively stronger statements: the gap shrinks, the gap is shrinking
% at most instants, and the CURRENT velocities already point at a near miss.
range=hypot(lure(:,1)-torso(:,1),lure(:,2)-torso(:,2));
closingSpeed=-trailingVelocity(range,t,fs,opts.velocityWindowSec);  % +ve = gap shrinking
% Prospective miss: distance at closest approach if both held their current
% velocity from here on.  Purely causal - trailing-window velocities only.
% This one is a statement about SUSTAINED motion, so it uses the longer
% missVelocityWindowSec: at the 100 ms window the relative-velocity estimate is
% noisy enough that tStar flickers about zero and the miss distance jumps
% between its two branches, which is estimator chatter rather than behaviour.
wMiss=opts.missVelocityWindowSec;
vLureX=trailingVelocity(lure(:,1),t,fs,wMiss);
vLureY=trailingVelocity(lure(:,2),t,fs,wMiss);
vCx=trailingVelocity(torso(:,1),t,fs,wMiss);
vCy=trailingVelocity(torso(:,2),t,fs,wMiss);
relX=lure(:,1)-torso(:,1);relY=lure(:,2)-torso(:,2);
relVX=vLureX-vCx;relVY=vLureY-vCy;
speed2=relVX.^2+relVY.^2;
tStar=-(relX.*relVX+relY.*relVY)./max(speed2,eps);
tStar(~isfinite(tStar))=NaN;
missDistance=range;                        % already past closest approach
ahead=isfinite(tStar)&tStar>0;
missDistance(ahead)=hypot(relX(ahead)+relVX(ahead).*tStar(ahead), ...
                          relY(ahead)+relVY(ahead).*tStar(ahead));

gRange=valid&isfinite(range);
firstIx=find(gRange,1,'first');lastIx=find(gRange,1,'last');
m.rangeStartM =range(firstIx);
m.rangeEndM   =range(lastIx);
m.rangeChangeM=range(lastIx)-range(firstIx);          % <0 = closed on the lure
gClose=valid&isfinite(closingSpeed);
m.medClosingSpeedMps=median(closingSpeed(gClose));    % >0 = closing
m.fracFramesClosing =mean(closingSpeed(gClose)>0);
gMiss=valid&isfinite(missDistance);
m.medMissDistanceM=median(missDistance(gMiss));
% Trend as last-third minus first-third medians, not endpoint minus endpoint:
% the endpoints of a jagged series are the least reliable samples in it.
m.missDistanceChangeM=thirdsTrend(missDistance,gMiss);   % <0 = aim improving
m.medTimeToClosestSec=median(tStar(gMiss&ahead),'omitnan');
m.fracFramesApproaching=mean(ahead(gMiss));   % current velocities still converging

% ---- directional coupling, strictly disjoint predictor/target windows ------
nTau=numel(opts.tauGridSec);
m.pcHeadTurnToTorso    =nan(1,nTau);
m.pcTorsoTurnToHead    =nan(1,nTau);
m.pcEtaToTorsoTurn     =nan(1,nTau);
m.pcHeadLeadToTorsoTurn=nan(1,nTau);
for q=1:nTau
    n=max(1,round(opts.tauGridSec(q)*fs));
    idx=(1+n):(N-n);
    if numel(idx)<opts.minPairs,continue;end
    ok=valid(idx)&valid(idx+n)&valid(idx-n);

    turnTorsoNext=uTorso(idx+n)-uTorso(idx);   % target  window [t, t+tau]
    turnTorsoPrev=uTorso(idx)-uTorso(idx-n);   % own momentum   [t-tau, t]
    turnHeadPrev =uHead(idx)-uHead(idx-n);     % predictor      [t-tau, t]
    turnHeadNext =uHead(idx+n)-uHead(idx);
    dLamTorso    =uLamTorso(idx+n)-uLamTorso(idx);  % common drive over target window
    dLamHead     =uLamHead(idx+n)-uLamHead(idx);
    turnHeadPrevOwn=turnHeadPrev;
    etaDev  =deviation(eta(idx));
    leadDev =deviation(deltaHead(idx));
    dTorsoDev=deviation(deltaTorso(idx));

    % head turn -> torso turn, given torso momentum and the LOS it is chasing
    m.pcHeadTurnToTorso(q)=partialCorr(turnHeadPrev,turnTorsoNext, ...
        [turnTorsoPrev,dLamTorso],ok,opts.minPairs);
    % reverse: torso turn -> head turn, symmetric conditioning
    m.pcTorsoTurnToHead(q)=partialCorr(turnTorsoPrev,turnHeadNext, ...
        [turnHeadPrevOwn,dLamHead],ok,opts.minPairs);
    % gaze offset deviation -> torso catch-up turn
    m.pcEtaToTorsoTurn(q)=partialCorr(etaDev,turnTorsoNext, ...
        [turnTorsoPrev,dTorsoDev],ok,opts.minPairs);
    % head's own lead deviation -> torso catch-up turn
    m.pcHeadLeadToTorsoTurn(q)=partialCorr(leadDev,turnTorsoNext, ...
        [turnTorsoPrev,dTorsoDev],ok,opts.minPairs);
end

% ---- effective sample size of the torso turn series ------------------------
n1=max(1,round(opts.tauGridSec(1)*fs));
m.effectiveNTorso=effectiveN(uTorso(1+n1:end)-uTorso(1:end-n1));

% ---- per-frame nominal series for plotting --------------------------------
m.frames=Drun.frames;m.t=t;m.valid=valid;
m.psiHeadDeg=rad2deg(psiHead);m.psiTorsoDeg=rad2deg(psiTorso);
m.lamHeadDeg=rad2deg(lamHead);m.lamTorsoDeg=rad2deg(lamTorso);
m.travelDeg=rad2deg(travel);
m.etaDeg=rad2deg(eta);m.deltaHeadDeg=rad2deg(deltaHead);
m.deltaTorsoDeg=rad2deg(deltaTorso);
m.rangeTorsoM=range;
m.travelMinusTorsoDeg=rad2deg(travelMinusTorso);
m.travelMinusHeadDeg =rad2deg(travelMinusHead);
m.closingSpeedMps=closingSpeed;
m.missDistanceM=missDistance;
end

% =========================================================================
function T=makeSurrogateTable(D,nominal,names,labels,opts)
% Circular-shift null for the head->torso partial correlation.  Shifting the head
% series relative to the torso series preserves each one's autocorrelation and
% destroys only their alignment, so the resulting |r| distribution is what a run
% this short and this smooth produces with no coupling whatsoever.
rows={};
for j=1:numel(names)
    d=D(j);pose=d.nominalPose;lure=d.nominalLure;fs=d.fs;N=size(pose,1);
    psiHead =forwardYaw(pose(:,1),pose(:,2),pose(:,3));
    psiTorso=forwardYaw(pose(:,7),pose(:,8),pose(:,9));
    lamTorso=atan2(lure(:,2)-pose(:,11),lure(:,1)-pose(:,10));
    uHead=unwrap(psiHead);uTorso=unwrap(psiTorso);uLamTorso=unwrap(lamTorso);
    valid=d.analysisMask(:);
    for q=1:numel(opts.tauGridSec)
        n=max(1,round(opts.tauGridSec(q)*fs));
        idx=(1+n):(N-n);
        if numel(idx)<opts.minPairs,continue;end
        ok=valid(idx)&valid(idx+n)&valid(idx-n);
        turnTorsoNext=uTorso(idx+n)-uTorso(idx);
        turnTorsoPrev=uTorso(idx)-uTorso(idx-n);
        dLamTorso=uLamTorso(idx+n)-uLamTorso(idx);
        headSeries=uHead(idx)-uHead(idx-n);
        obs=partialCorr(headSeries,turnTorsoNext,[turnTorsoPrev,dLamTorso], ...
            ok,opts.minPairs);
        nShift=min(opts.nSurrogate,numel(idx)-1);
        shifts=round(linspace(1,numel(idx)-1,nShift));
        null=nan(numel(shifts),1);
        for s=1:numel(shifts)
            null(s)=partialCorr(circshift(headSeries,shifts(s)),turnTorsoNext, ...
                [turnTorsoPrev,dLamTorso],ok,opts.minPairs);
        end
        pTwoSided=mean(abs(null)>=abs(obs),'omitnan');
        rows{end+1}=table(string(names{j}),string(labels{j}), ...
            1000*opts.tauGridSec(q),obs,percentile(abs(null),95), ...
            median(abs(null),'omitnan'),pTwoSided,nominal(j).effectiveNTorso, ...
            'VariableNames',{'run','label','tauMs','observedPartialR', ...
            'null_absR_p95','null_absR_median','surrogateP','effectiveN'}); %#ok<AGROW>
    end
end
T=vertcat(rows{:});
end

% =========================================================================
function T=makeRunSummary(nominal,mc,names,labels,animals,flds)
nRun=numel(names);rows=cell(nRun,1);
for j=1:nRun
    row=table(string(names{j}),string(labels{j}),string(animals{j}), ...
        nominal(j).N,'VariableNames',{'run','label','animal','nFrames'});
    for f=1:numel(flds)
        col=mc.(flds{f})(:,j);
        row.(['nom_' flds{f}])=nominal(j).(flds{f});
        row.(['mc_' flds{f} '_median'])=percentile(col,50);
        row.(['mc_' flds{f} '_lo'])=percentile(col,2.5);
        row.(['mc_' flds{f} '_hi'])=percentile(col,97.5);
    end
    row.pTorsoLockedTighterThanHead=mean(mc.lockRatioTorsoOverHead(:,j)<1);
    rows{j}=row;
end
T=vertcat(rows{:});
end

function T=makeCouplingTable(nominal,mc,names,labels,animals,vecFlds,opts)
rows={};
for j=1:numel(names)
    for q=1:numel(opts.tauGridSec)
        row=table(string(names{j}),string(labels{j}),string(animals{j}), ...
            1000*opts.tauGridSec(q),'VariableNames',{'run','label','animal','tauMs'});
        for f=1:numel(vecFlds)
            col=mc.(vecFlds{f})(:,j,q);
            row.(['nom_' vecFlds{f}])=nominal(j).(vecFlds{f})(q);
            row.(['mc_' vecFlds{f} '_median'])=percentile(col,50);
            row.(['mc_' vecFlds{f} '_lo'])=percentile(col,2.5);
            row.(['mc_' vecFlds{f} '_hi'])=percentile(col,97.5);
            row.(['pPositive_' vecFlds{f}])=mean(col>0);
        end
        rows{end+1}=row; %#ok<AGROW>
    end
end
T=vertcat(rows{:});
end

function T=makeBetweenRunSigns(mc,vecFlds,opts)
% The honest verdict layer: n=5 runs, not n=frames.  For a scalar metric the sign
% test is on the metric itself; for the coupling correlations it is on the
% per-run MC median at each tau.  fracDrawsAllPositive additionally asks how
% often a single joint reconstruction makes all five runs agree.
rows={};
% Ratio metrics are tested against 1 on the log scale, oriented so that
% "positive" always means "the effect claimed by the hypothesis".
rows{end+1}=signRow('lockRatioTorsoOverHead_lt1',NaN, ...
    -log(mc.lockRatioTorsoOverHead));
rows{end+1}=signRow('losLockGainTorso_gt1',NaN,log(mc.losLockGainTorso));
rows{end+1}=signRow('losLockGainHead_gt1',NaN,log(mc.losLockGainHead));
% Oriented so "positive" = the claim: torso closer to travel than head is; the
% gap to the lure shrinks; the animal is closing; the prospective miss improves.
rows{end+1}=signRow('travelAlignmentGainDeg_gt0',NaN,mc.travelAlignmentGainDeg);
rows{end+1}=signRow('rangeChange_lt0',NaN,-mc.rangeChangeM);
rows{end+1}=signRow('medClosingSpeed_gt0',NaN,mc.medClosingSpeedMps);
rows{end+1}=signRow('fracFramesClosing_gt_half',NaN,mc.fracFramesClosing-0.5);
rows{end+1}=signRow('missDistanceChange_lt0',NaN,-mc.missDistanceChangeM);
for i=1:numel(vecFlds)
    for q=1:numel(opts.tauGridSec)
        X=squeeze(mc.(vecFlds{i})(:,:,q));
        rows{end+1}=signRow(vecFlds{i},1000*opts.tauGridSec(q),X); %#ok<AGROW>
    end
end
T=vertcat(rows{:});
end

function row=signRow(metric,tauMs,X)
med=median(X,1,'omitnan');med=med(isfinite(med));
nPos=sum(med>0);nRun=numel(med);
row=table(string(metric),tauMs,nPos,nRun,signTestP(nPos,nRun), ...
    mean(all(X>0,2)),mean(all(X<0,2)),mean(med), ...
    'VariableNames',{'metric','tauMs','nRunsPositive','nRuns','signTestP', ...
    'fracDrawsAllPositive','fracDrawsAllNegative','meanAcrossRuns'});
end

function T=makeByFrameTable(nominal,names,animals)
rows=cell(numel(nominal),1);
for j=1:numel(nominal)
    m=nominal(j);n=numel(m.frames);
    rows{j}=table(repmat(string(names{j}),n,1),repmat(string(animals{j}),n,1), ...
        m.frames(:),m.t(:),double(m.valid(:)),m.psiHeadDeg(:),m.psiTorsoDeg(:), ...
        m.lamHeadDeg(:),m.lamTorsoDeg(:),m.travelDeg(:),m.etaDeg(:), ...
        m.deltaHeadDeg(:),m.deltaTorsoDeg(:),m.rangeTorsoM(:), ...
        m.travelMinusHeadDeg(:),m.travelMinusTorsoDeg(:), ...
        m.closingSpeedMps(:),m.missDistanceM(:), ...
        'VariableNames',{'run','animal','absFrame','timeSec','valid', ...
        'headYawDeg','torsoYawDeg','headLosBearingDeg','torsoLosBearingDeg', ...
        'travelHeadingDeg','gazeOffsetDeg','headLeadDeg','torsoOffsetDeg', ...
        'torsoRangeM','travelMinusHeadDeg','travelMinusTorsoDeg', ...
        'closingSpeedMps','missDistanceM'});
end
T=vertcat(rows{:});
end

% =========================================================================
function printHeadline(runSummary,couplingTau,betweenRun,opts)
fprintf('\n---- premise: is the separation real? ----\n');
for j=1:height(runSummary)
    fprintf('%-18s |travel-torso| %4.1f | |gaze offset| %4.1f | head-to-lure %4.1f | torso-to-lure %4.1f deg\n', ...
        runSummary.label(j),runSummary.nom_medAbsTravelMinusTorsoDeg(j), ...
        runSummary.nom_medAbsEtaDeg(j),runSummary.nom_medAbsDeltaHeadDeg(j), ...
        runSummary.nom_medAbsDeltaTorsoDeg(j));
end
fprintf('\n---- which axis is the locomotor axis? (deg from travel heading) ----\n');
for j=1:height(runSummary)
    fprintf('%-18s head %5.1f | torso %5.1f | torso closer by %5.1f deg\n', ...
        runSummary.label(j),runSummary.nom_medAbsTravelMinusHeadDeg(j), ...
        runSummary.nom_medAbsTravelMinusTorsoDeg(j), ...
        runSummary.nom_travelAlignmentGainDeg(j));
end
fprintf('\n---- does the locomotion close on the lure? ----\n');
for j=1:height(runSummary)
    fprintf('%-18s range %4.2f -> %4.2f m (%+5.2f) | closing %+5.2f m/s, %3.0f%% of frames | prospective miss %4.2f m (%+5.2f)\n', ...
        runSummary.label(j),runSummary.nom_rangeStartM(j), ...
        runSummary.nom_rangeEndM(j),runSummary.nom_rangeChangeM(j), ...
        runSummary.nom_medClosingSpeedMps(j),100*runSummary.nom_fracFramesClosing(j), ...
        runSummary.nom_medMissDistanceM(j),runSummary.nom_missDistanceChangeM(j));
end
fprintf('\n---- LOS lock: offset +/- RMS scatter, vs trivial nulls (deg) ----\n');
fprintf('   lock gain >1 means the LOS-referenced offset beats holding a constant world heading\n');
for j=1:height(runSummary)
    fprintf('%-18s head %5.1f +/-%4.1f (gain %4.2f) | torso %5.1f +/-%4.1f (gain %4.2f) | ratio %4.2f, P<1=%.2f\n', ...
        runSummary.label(j),runSummary.nom_deltaHeadMeanDeg(j), ...
        runSummary.nom_deltaHeadRmsDeg(j),runSummary.nom_losLockGainHead(j), ...
        runSummary.nom_deltaTorsoMeanDeg(j),runSummary.nom_deltaTorsoRmsDeg(j), ...
        runSummary.nom_losLockGainTorso(j), ...
        runSummary.nom_lockRatioTorsoOverHead(j), ...
        runSummary.pTorsoLockedTighterThanHead(j));
end
fprintf('\n---- power ceiling: autocorrelation-effective sample size per run ----\n');
for j=1:height(runSummary)
    fprintf('%-18s %d frames -> effective n = %.1f\n',runSummary.label(j), ...
        runSummary.nFrames(j),runSummary.nom_effectiveNTorso(j));
end
fprintf('\n---- coupling (MC median partial r) by run x tau ----\n');
for q=1:numel(opts.tauGridSec)
    tt=1000*opts.tauGridSec(q);
    sub=couplingTau(couplingTau.tauMs==tt,:);
    fprintf('tau=%3.0f ms  head->torso: ',tt);
    fprintf('%6.2f',sub.mc_pcHeadTurnToTorso_median);
    fprintf('   |  torso->head: ');
    fprintf('%6.2f',sub.mc_pcTorsoTurnToHead_median);
    fprintf('\n');
end
fprintf('\n---- between-run sign tests (n=5 runs) ----\n');
disp(betweenRun);
end

% =========================================================================
function y=forwardYaw(tx,ty,tz)
% World yaw of the segment-forward axis R*[1;0;0] with R=Rz(tz)*Rx(tx)*Ry(ty)
% under the CRBR rotation convention (identical algebra for head and torso).
hx=cos(tz).*cos(ty)-sin(tz).*sin(tx).*sin(ty);
hy=-sin(tz).*cos(ty)-cos(tz).*sin(tx).*sin(ty);
y=atan2(hy,hx);
end

function v=trailingVelocity(x,t,fs,windowSec)
% Causal trailing quadratic fit, velocity at the window end.  Used only for the
% travel-heading premise check, never as a predictor.
n=numel(x);v=nan(n,1);w=max(4,round(windowSec*fs));
for i=w:n
    ix=(i-w+1):i;
    if any(~isfinite(x(ix))),continue;end
    tt=t(ix)-t(i);X=[ones(w,1),tt(:),tt(:).^2];b=X\x(ix);v(i)=b(2);
end
end

function d=thirdsTrend(x,good)
% Last-third median minus first-third median, over the valid samples only.
% Robust to the spikes that endpoint-minus-endpoint would be dominated by.
ix=find(good(:)&isfinite(x(:)));
if numel(ix)<6,d=NaN;return;end
k=floor(numel(ix)/3);
d=median(x(ix(end-k+1:end)))-median(x(ix(1:k)));
end

function y=deviation(a)
% Signed departure from the run's own circular mean.
b=circMean(a);y=wrapPi(a-b);
end

function r=partialCorr(x,y,ctrl,ok,minPairs)
x=x(:);y=y(:);
g=ok(:)&isfinite(x)&isfinite(y)&all(isfinite(ctrl),2);
if nnz(g)<minPairs,r=NaN;return;end
Z=[ones(nnz(g),1),ctrl(g,:)];
if rank(Z)<size(Z,2),r=NaN;return;end
rx=x(g)-Z*(Z\x(g));ry=y(g)-Z*(Z\y(g));
if std(rx)<1e-12||std(ry)<1e-12,r=NaN;else,r=corr(rx,ry);end
end

function n=effectiveN(x)
% Autocorrelation-corrected sample size, n*(1-r1)/(1+r1) with lag-1 r1.
x=x(isfinite(x));n=numel(x);
if n<4,return;end
x=x-mean(x);r1=sum(x(1:end-1).*x(2:end))/max(sum(x.^2),eps);
r1=min(max(r1,-0.99),0.99);n=n*(1-r1)/(1+r1);
end

function c=circMean(a),a=a(isfinite(a));if isempty(a),c=NaN;else,c=atan2(mean(sin(a)),mean(cos(a)));end,end
function v=rmsAboutCircMean(a)
a=a(isfinite(a));if isempty(a),v=NaN;return;end
v=rad2deg(sqrt(mean(wrapPi(a-circMean(a)).^2)));
end
function p=signTestP(k,n)
if n<1,p=NaN;return;end
kk=max(k,n-k);p=0;for i=kk:n,p=p+nchoosek(n,i)*0.5^n;end,p=min(1,2*p);
end
function p=percentile(x,q)
x=sort(x(isfinite(x)));if isempty(x),p=NaN;return;end
u=1+(numel(x)-1)*q/100;lo=floor(u);hi=ceil(u);
if lo==hi,p=x(lo);else,p=x(lo)+(u-lo)*(x(hi)-x(lo));end
end
function x=wrapPi(x),x=mod(x+pi,2*pi)-pi;end
function opts=applyDefaults(opts,defaults)
names=fieldnames(defaults);
for i=1:numel(names)
    if ~isfield(opts,names{i})||isempty(opts.(names{i}))
        opts.(names{i})=defaults.(names{i});
    end
end
end
