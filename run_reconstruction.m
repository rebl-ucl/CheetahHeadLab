function result = run_reconstruction(dataset, varargin)
%RUN_RECONSTRUCTION Reconstruct one pursuit from released image annotations.
% run_reconstruction()                         % complete 69-frame example
% run_reconstruction('julesTopRun1_20171209')   % another released pursuit
% run_reconstruction(name,'nMC',500)           % also reconstruct uncertainty
% run_reconstruction(name,'inputFolder',path)  % self-calibration output
% All work is written to outputRoot/dataset. Inputs and references are read-only.

root = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(dataset), dataset = 'julesBottomFlick2_20171209'; end
dataset = char(dataset);
p = inputParser;
p.addParameter('inputFolder','',@(x)ischar(x)||isstring(x));
p.addParameter('outputRoot',fullfile(root,'outputs'),@(x)ischar(x)||isstring(x));
p.addParameter('nMC',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x==fix(x));
p.addParameter('useParallel',false,@islogical);
p.parse(varargin{:}); opts=p.Results;
M = jsondecode(fileread(fullfile(root,'manifest.json')));
idx = find(strcmp({M.dataset},dataset));
assert(isscalar(idx),'Choose a dataset listed in manifest.json.');
entry=M(idx);
inputFolder=char(opts.inputFolder);
if isempty(inputFolder), inputFolder=fullfile(root,'data','inputs',dataset); end
work=fullfile(char(opts.outputRoot),dataset);
assert(isfolder(inputFolder),'Input folder not found: %s',inputFolder);
assert(~isfolder(work),'Output already exists: %s. Choose a new outputRoot to rerun.',work);

oldPath=path;
oldVisibility=get(groot,'defaultFigureVisible');
cleanup=onCleanup(@()restoreSession(oldPath,oldVisibility)); %#ok<NASGU>
addpath(fullfile(root,'CRBR','reconstruction'), ...
    fullfile(root,'CRBR','uncertainty'),fullfile(root,'CRBR','lure'));
set(groot,'defaultFigureVisible','off');
if ~isfolder(opts.outputRoot), mkdir(opts.outputRoot); end
copyfile(inputFolder,work);
diary(fullfile(work,'reconstruction.log'));
logCleanup=onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('CRBR reconstruction: %s, %g frames/s\n',dataset,entry.fps);
t=tic;
step=struct('showFigures',false,'saveFigures',false,'showVideos',false, ...
    'saveVideos',false,'createVideo',false,'createOrthoVideo',false, ...
    'plotCameraConfig',false,'fs',entry.fps,'optimSuffix',"", ...
    'robustOpt',struct('type','huber','c',3));
camIDs=entry.camera_ids(:)';
multiviewRecon(work,camIDs,step);
estimateScales(work,step);
reconPerFrame(work,camIDs,step);
lureReconTriangulation(work,camIDs,step);
step.camIDs=camIDs;
step.outSuffix="_trajectory";
result=reconTrajectoryBatch(work,step);
assert(all(isfinite(result.pose),'all'),'Reconstruction contains non-finite pose parameters.');

% Export readable coordinates without requiring MATLAB to inspect the output.
D=mcLoadDataset(work,camIDs);
names={'head_rx_rad','head_ry_rad','head_rz_rad','head_x_m','head_y_m','head_z_m', ...
    'torso_rx_rad','torso_ry_rad','torso_rz_rad','torso_x_m','torso_y_m','torso_z_m'};
T=array2table(result.pose,'VariableNames',names);
T=addvars(T,result.absFrames,'Before',1,'NewVariableNames','absolute_frame');
writetable(T,fullfile(work,'poses.csv'));
Rx=@(a)[1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)];
Ry=@(a)[cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)];
Rz=@(a)[cos(a) sin(a) 0;-sin(a) cos(a) 0;0 0 1];
ang=nan(numel(result.absFrames),6);
for k=1:size(ang,1)
    ang(k,:)=rad2deg(mcMainAngles(result.pose(k,:),D.lureXYZ(k,:),Rz,Rx,Ry));
end
A=array2table([result.absFrames D.lureXYZ ang], 'VariableNames', ...
    {'absolute_frame','lure_x_m','lure_y_m','lure_z_m','alpha_head_pitch_error_deg', ...
    'beta_torso_pitch_error_deg','gamma_head_yaw_error_deg','delta_torso_yaw_error_deg', ...
    'head_lure_angle_deg','torso_lure_angle_deg'});
writetable(A,fullfile(work,'angles_and_lure.csv'));
summary=struct('dataset',dataset,'fps',entry.fps,'frames',numel(result.absFrames), ...
    'elapsed_seconds',toc(t),'matlab_version',version,'input_folder',inputFolder);

% Report forward-axis and origin differences; Euler triples can be ambiguous.
ref=load(fullfile(root,'reference',dataset,'storedOptimRaw_trajectory.mat'));
[common,ia,ib]=intersect(result.absFrames,ref.absFrames);
summary.reference_shared_frames=numel(common);
summary.reference_frames=numel(ref.absFrames);
summary.identical_frame_coverage=isequal(result.absFrames(:),ref.absFrames(:));
if ~isempty(common)
    cur=result.pose(ia,:); old=ref.xBarStore(ib,:);
    err=nan(numel(common),4);
    for k=1:numel(common)
        for body=0:1
            q=body*6;
            h=Rz(cur(k,q+3))*Rx(cur(k,q+1))*Ry(cur(k,q+2))*[1;0;0];
            hr=Rz(old(k,q+3))*Rx(old(k,q+1))*Ry(old(k,q+2))*[1;0;0];
            err(k,body+1)=atan2d(norm(cross(h,hr)),dot(h,hr));
            err(k,body+3)=norm(cur(k,q+(4:6))-old(k,q+(4:6)));
        end
    end
    C=array2table([common err],'VariableNames',{'absolute_frame','head_axis_difference_deg', ...
        'torso_axis_difference_deg','head_origin_difference_m','torso_origin_difference_m'});
    writetable(C,fullfile(work,'reference_comparison.csv'));
    summary.max_head_axis_difference_deg=max(err(:,1));
    summary.max_torso_axis_difference_deg=max(err(:,2));
    summary.max_head_origin_difference_m=max(err(:,3));
    summary.max_torso_origin_difference_m=max(err(:,4));
end
fid=fopen(fullfile(work,'reproduction_summary.json'),'w');
assert(fid>=0,'Could not write reproduction summary.');
fprintf(fid,'%s\n',jsonencode(summary,PrettyPrint=true));fclose(fid);
disp(summary);
if opts.nMC>0
    run_monte_carlo(work,opts.nMC,opts.useParallel);
end
fprintf('Finished. Outputs: %s\n',work);
end

function restoreSession(oldPath,oldVisibility)
path(oldPath);
set(groot,'defaultFigureVisible',oldVisibility);
end
