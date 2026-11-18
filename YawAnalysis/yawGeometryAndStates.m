function yawGeometryAndStates(C,out)
%YAWGEOMETRYANDSTATES Propagate joint reconstructions into geometry and states.
% The equal-range projection sets ahead/lagged direction, independently of
% which side of the cheetah the lure occupies. Aligned draws use |gamma|<=gate.
if ~isfolder(out), mkdir(out); end
cfg=struct('pathWindowSec',.10,'pathMaxWindowSec',.30,'pathMinDisplacementM',.05, ...
    'pathMinLinearity',.65,'hampelK',3,'hampelSigma',3);
parts={}; runs={}; draws={}; gates={};
contexts=["pitch tracking","steady lateral","close pass-by"];
contextNames=["frontal pursuit","lateral pursuit","close pass-by pursuit"];
for j=1:height(C.configuration)
    row=C.configuration(j,:); folder=char(row.dataDir); name=char(row.dataset);
    A=load(fullfile(folder,'lureAnalysis_results_trajectory.mat'),'runResults');
    S=load(fullfile(folder,'MonteCarloResults',['mcMain_' name '.mat']), ...
        'absFrames','poseAll','lureXYZAll','sampleSucceeded');
    frames=double(S.absFrames(:)); n=numel(frames); fs=A.runResults.fs;
    assert(isequal(frames,C.frames.absFrame(C.frames.dataset==row.dataset)));
    assert(all(diff(frames)==1)); t=(frames-frames(1))/fs;
    ok=find(S.sampleSucceeded); nd=numel(ok); assert(nd==500);
    gamma=nan(nd,n); lead=gamma; theta=gamma; range=gamma; gain=gamma; borrowed=false(nd,n);
    for m=1:nd
        Q=yawSampleGeometry(double(squeeze(S.poseAll(ok(m),:,:))), ...
            double(squeeze(S.lureXYZAll(ok(m),:,:))),t,true(n,1),fs,cfg);
        r=Q.lureXYZ(:,1:2)-Q.eyeXYZ(:,1:2); distance=vecnorm(r,2,2);
        h=[cos(Q.headYawRad),sin(Q.headYawRad)];
        score=sum((h-r./distance).*Q.pathTangent,2);
        assert(all(Q.valid)&&all(distance>0));
        assert(max(abs(score-Q.pathLeadM./distance))<1e-10);
        theta(m,:)=abs(rad2deg(mod(atan2(Q.pathTangent(:,2),Q.pathTangent(:,1))-Q.bearingCurrentRad+pi,2*pi)-pi));
        gamma(m,:)=Q.gammaDeg; lead(m,:)=Q.pathLeadM; range(m,:)=distance;
        gain(m,:)=score; borrowed(m,:)=Q.directionBorrowed;
    end
    signed=sign(lead).*abs(gamma);
    point=sign(median(lead,1).').*abs(median(gamma,1).');
    thm=median(theta,1).'; rm=median(range,1).';
    rt=median(theta,2); rr=median(range,2); ry=median(signed,2);
    qt=prctile(rt,[2.5 97.5]); qr=prctile(rr,[2.5 97.5]); qy=prctile(ry,[2.5 97.5]);
    context=contexts(contextNames==row.context); [group,~]=fileparts(C.datasets{j});
    parts{j}=table(repmat(row.dataset,n,1),repmat(string(group),n,1),repmat(context,n,1), ...
        frames,point,thm,rm,mean(borrowed,1).','VariableNames', ...
        {'dataset','group','context','absFrame','signedPathGammaDeg','thetaMedianDeg','rangeMedianM','directionBorrowedFraction'});
    runs{j}=table(row.dataset,string(group),context,lower(row.animal),n,n,median(point),qy(1),qy(2), ...
        median(thm),qt(1),qt(2),median(rm),qr(1),qr(2),prctile(thm,25),prctile(thm,75), ...
        prctile(rm,25),prctile(rm,75),mean(borrowed,'all'),fs,nd,'VariableNames', ...
        {'dataset','group','context','animal','nYawFrames','nGeometryFrames','medianSignedYawDeg', ...
        'yawReconLoDeg','yawReconHiDeg','medianThetaDeg','thetaReconLoDeg','thetaReconHiDeg', ...
        'medianRangeM','rangeReconLoM','rangeReconHiM','thetaFrameQ25Deg','thetaFrameQ75Deg', ...
        'rangeFrameQ25M','rangeFrameQ75M','directionBorrowedFraction','fsHz','nDraws'});
    draws{j}=table(repmat(row.dataset,nd,1),(1:nd).',ok,rt,rr,ry,'VariableNames', ...
        {'dataset','draw','successfulSampleIndex','medianThetaDeg','medianRangeM','medianSignedYawDeg'});
    for gate=[5 10 15]
        aligned=abs(gamma)<=gate; ahead=~aligned&gain>1e-10; lagged=~aligned&gain< -1e-10;
        ambiguous=~aligned&abs(gain)<=1e-10;
        assert(all(double(aligned)+ahead+lagged+ambiguous==1,'all'));
        pA=mean(ahead,1).'; pZ=mean(aligned,1).'; pL=mean(lagged,1).'; pU=mean(ambiguous,1).';
        state=repmat("unresolved",n,1); state(pA>=.95)="ahead"; state(pZ>=.95)="aligned"; state(pL>=.95)="lagged";
        gates{end+1}=table(repmat(row.dataset,n,1),repmat(row.animal,n,1),repmat(row.context,n,1), ...
            frames,repmat(gate,n,1),pA,pZ,pL,pU,state,'VariableNames', ...
            {'dataset','animal','context','absFrame','gateDeg','pAhead','pAligned','pLagged','pDirectionAmbiguous','supportedState'}); %#ok<AGROW>
    end
    fprintf('Yaw geometry: %s, %d frames, %d draws\n',name,n,nd);
end
G=vertcat(parts{:}); R=vertcat(runs{:}); D=vertcat(draws{:}); T=vertcat(gates{:});
assert(height(G)==1708 && height(R)==18 && height(T)==3*1708);
writetable(G,fullfile(out,'geometry_by_frame.csv'));
writetable(G(:,1:5),fullfile(out,'signed_yaw_by_frame_all18.csv'));
writetable(R,fullfile(out,'geometry_by_run.csv')); writetable(D,fullfile(out,'geometry_run_draws.csv'));
writetable(T,fullfile(out,'gate_support_by_frame.csv'));
end
