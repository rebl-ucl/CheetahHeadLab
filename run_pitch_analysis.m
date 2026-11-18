function run_pitch_analysis(sourceRoot)
%RUN_PITCH_ANALYSIS Manuscript pitch statistics for 16 pursuits / 1590 frames.
% Run reconstruction for each pitch-included pursuit first.
if nargin<1, sourceRoot=[]; end
C=prepare_analysis_inputs(sourceRoot,'pitch',false);
old=path; cleanup=onCleanup(@()path(old)); %#ok<NASGU>
addpath(fullfile(C.root,'PitchAnalysis'));
base=fullfile(C.root,'outputs','pitch');
opts=struct('dataRoot',C.dataRoot,'showFigures',false,'saveFigures',false, ...
    'useDespiked',false,'saveCounterRotationOutputs',false,'showFigureTitles',false);
opts.outputDir=fullfile(base,'tightness_run_equal'); opts.ecdfWeighting='run';
pitchTightnessAnalysis(C.datasets,opts);
opts.outputDir=fullfile(base,'tightness_animal_equal'); opts.ecdfWeighting='animal'; opts.saveFigures=true;
pitchTightnessAnalysis(C.datasets,opts);
pitchControlModelAnalysis(C.datasets,struct('dataRoot',C.dataRoot, ...
    'outputDir',fullfile(base,'mechanism'),'showFigures',false,'saveFigures',false));
writetable(C.configuration,fullfile(base,'cohort.csv'));
fprintf('Pitch analysis: 16 pursuits, 1590 frames. Figures/tables: %s\n',base);
end
