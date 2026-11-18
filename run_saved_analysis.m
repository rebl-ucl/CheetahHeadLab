function run_saved_analysis(secondary)
%RUN_SAVED_ANALYSIS Recompute manuscript statistics from bundled joint draws.
% run_saved_analysis        % includes supplementary spatial/closure analyses
% run_saved_analysis(false) % main pitch/yaw analyses only
if nargin<1, secondary=true; end
root=fileparts(mfilename('fullpath'));
source=fullfile(root,'data','analysis_cache');
assert(isfolder(source),'Missing bundled data/analysis_cache.');
run_pitch_analysis(source);
run_yaw_analysis(source,secondary);
fprintf('Next: python PitchAnalysis/plot_pitch.py\n');
fprintf('      python YawAnalysis/plot_yaw.py\n');
fprintf('      python YawAnalysis/compare_models.py\n');
end
