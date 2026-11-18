function result = run_monte_carlo(workFolder,nDraws,useParallel)
%RUN_MONTE_CARLO Refit joint uncertainty draws after run_reconstruction.
% run_monte_carlo('outputs/julesBottomFlick2_20171209',500)
% Two draws can check execution, but cannot estimate uncertainty intervals.
if nargin<2, nDraws=500; end
if nargin<3, useParallel=false; end
validateattributes(nDraws,{'numeric'},{'scalar','integer','positive'});
root=fileparts(mfilename('fullpath'));
oldPath=path;cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'CRBR','uncertainty'),fullfile(root,'CRBR','reconstruction'));
assert(isfile(fullfile(workFolder,'storedOptimRaw_trajectory.mat')), ...
    'Run reconstruction first.');
[ok,info]=fileattrib(workFolder);assert(ok,'Output folder not found.');
[parent,dataset]=fileparts(info.Name);
M=jsondecode(fileread(fullfile(root,'manifest.json')));
entry=M(strcmp({M.dataset},dataset));
assert(isscalar(entry),'Unknown dataset.');
out=fullfile(info.Name,sprintf('MonteCarlo_%d',nDraws));
assert(~isfolder(out),'MC output exists: %s. Use a fresh output folder to rerun.',out);
opts=struct('dataRoot',parent,'camIDs',entry.camera_ids(:)', ...
    'nMC',nDraws,'seed',1,'makePlots',false,'useParallel',useParallel, ...
    'reconstructionMode','trajectory','outDir',out);
result=mcMainAnalysis(dataset,opts);
fprintf('Successful joint draws: %d / %d\n',sum(result.sampleSucceeded),nDraws);
if ~all(result.sampleSucceeded)
    warning('release:failedDraws', ...
        'Some reconstructions failed. Inspect sampleError before using uncertainty summaries.');
end
end
