function C = prepare_analysis_inputs(sourceRoot, cohort, requireMC)
%PREPARE_ANALYSIS_INPUTS Join CRBR outputs to the exact manuscript frame list.
% Sources are read only. Derived inputs go to outputs/analysis_inputs.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(sourceRoot), sourceRoot=fullfile(root,'outputs'); end
if nargin<2, cohort='yaw'; end
if nargin<3, requireMC=false; end
assert(ismember(string(cohort),["pitch","yaw"]));
M=jsondecode(fileread(fullfile(root,'manifest.json')));
F=readtable(fullfile(root,'reference','manuscript_frames.csv'),'TextType','string');
if strcmp(cohort,'pitch'), F=F(string(F.pitch_included)=="True" | string(F.pitch_included)=="1",:); end
viewRoot=fullfile(root,'outputs','analysis_inputs');
rows=cell(0,1); keys={};
for j=1:numel(M)
    e=M(j); name=e.dataset; selected=F.absolute_frame(F.dataset==string(name));
    if isempty(selected), continue; end
    folder=fullfile(sourceRoot,name);
    if ~isfolder(folder), folder=fullfile(sourceRoot,e.group,name); end
    poseFile=fullfile(folder,'storedOptimRaw_trajectory.mat');
    assert(isfile(poseFile),'Reconstruct %s first. Missing %s',name,poseFile);
    P=load(poseFile,'absFrames','xBarRawStore');
    L=load(fullfile(folder,'3DLurePos.mat'),'Lure3D');
    [ok,poseRows]=ismember(selected,double(P.absFrames(:))); assert(all(ok),'Missing selected poses: %s',name);
    [ok,li]=ismember(selected,double(L.Lure3D.Index(:))); assert(all(ok),'Missing selected lure positions: %s',name);
    p=double(P.xBarRawStore(poseRows,:)); marks=L.Lure3D.Mark;
    if iscell(marks), xyz=cell2mat(cellfun(@(x)reshape(double(x),1,[]),reshape(marks(li),[],1),'UniformOutput',false));
    else, xyz=double(marks(li,1:3)); end
    assert(all(isfinite(p),'all')&&all(isfinite(xyz),'all'));
    H=forward(p(:,1:3)); T=forward(p(:,7:9));
    [hy,hp]=cart2sph(H(:,1),H(:,2),H(:,3));
    [ty,tp]=cart2sph(T(:,1),T(:,2),T(:,3));
    dh=xyz-p(:,4:6); dt=xyz-p(:,10:12);
    [ly,lp,r]=cart2sph(dh(:,1),dh(:,2),dh(:,3));
    [by,bp,rb]=cart2sph(dt(:,1),dt(:,2),dt(:,3));
    wrap=@(x)mod(x+pi,2*pi)-pi;
    n=numel(selected);
    runResults=struct('commonFrames',selected,'nFrames',n,'fs',e.fps, ...
        'alpha',wrap(hp-lp),'beta',wrap(tp-bp),'gamma',wrap(hy-ly),'delta',wrap(ty-by), ...
        'theta_h_pitch',hp,'theta_t_pitch',tp,'lambda_h_pitch',lp,'lambda_t_pitch',bp, ...
        'theta_h_yaw',hy,'theta_t_yaw',ty,'lambda_h_yaw',ly,'lambda_t_yaw',by, ...
        'r',r,'rb',rb,'pitchAnalysisMask',true(n,1),'yawAnalysisMask',true(n,1), ...
        'isInterp',false(n,1),'despikeAngles',false);
    dest=fullfile(viewRoot,e.group,name); if ~isfolder(dest), mkdir(dest); end
    save(fullfile(dest,'lureAnalysis_results_trajectory.mat'),'runResults');
    copyfile(poseFile,fullfile(dest,'storedOptimRaw_trajectory.mat'));
    copyfile(fullfile(folder,'3DLurePos.mat'),fullfile(dest,'3DLurePos.mat'));
    pathFile=fullfile(folder,'3DLurePos_full.mat');
    if ~isfile(pathFile), pathFile=fullfile(folder,'3DLurePos.mat'); end
    copyfile(pathFile,fullfile(dest,'3DLurePos_full.mat'));
    if requireMC
        mcFile=fullfile(folder,'MonteCarlo_500',['mcMain_' name '.mat']);
        if ~isfile(mcFile), mcFile=fullfile(folder,'MonteCarloResults',['mcMain_' name '.mat']); end
        assert(isfile(mcFile),'Run 500 Monte Carlo draws first: %s',name);
        S=load(mcFile,'absFrames','poseAll','lureXYZAll','sampleSucceeded', ...
            'nominalPose','nominalLureXYZ','angleNames','medianDeg','opts','provenance');
        assert(nnz(S.sampleSucceeded)==500,'Need 500 successful joint draws: %s',name);
        [ok,mi]=ismember(selected,double(S.absFrames(:))); assert(all(ok));
        assert(max(abs(S.nominalPose(mi,:)-p),[],'all')<1e-8,'MC poses differ from nominal: %s',name);
        assert(max(abs(S.nominalLureXYZ(mi,:)-xyz),[],'all')<1e-8,'MC lure differs from nominal: %s',name);
        S.absFrames=selected; S.poseAll=S.poseAll(:,mi,:); S.lureXYZAll=S.lureXYZAll(:,mi,:);
        S.nominalPose=S.nominalPose(mi,:); S.nominalLureXYZ=S.nominalLureXYZ(mi,:); S.medianDeg=S.medianDeg(mi,:);
        mcDir=fullfile(dest,'MonteCarloResults'); if ~isfolder(mcDir), mkdir(mcDir); end
        save(fullfile(mcDir,['mcMain_' name '.mat']),'-struct','S','-v7.3');
    end
    context="frontal pursuit";
    if strcmp(e.group,'yaw_tracking'), context="lateral pursuit"; end
    if startsWith(name,'zorro'), context="close pass-by pursuit"; end
    animal=string(regexp(name,'^[a-z]+','match','once'));
    animal=upper(extractBefore(animal,2))+extractAfter(animal,1);
    rows{end+1}=table(string(name),animal,context,string(dest),n, ...
        'VariableNames',{'dataset','animal','context','dataDir','nFrames'}); %#ok<AGROW>
    keys{end+1,1}=[e.group '/' name]; %#ok<AGROW>
end
C=struct('root',root,'dataRoot',viewRoot,'datasets',{keys},'configuration',vertcat(rows{:}));
C.frames=table(F.dataset,F.absolute_frame,'VariableNames',{'dataset','absFrame'});
if strcmp(cohort,'pitch'), assert(height(C.configuration)==16 && height(F)==1590);
else, assert(height(C.configuration)==18 && height(F)==1708); end
end

function v=forward(a)
x=a(:,1); y=a(:,2); z=a(:,3);
v=[cos(z).*cos(y)-sin(z).*sin(x).*sin(y), ...
    -sin(z).*cos(y)-cos(z).*sin(x).*sin(y), -cos(x).*sin(y)];
end
