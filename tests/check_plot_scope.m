function check_plot_scope()
%CHECK_PLOT_SCOPE Direct calls cannot enable the removed analysis figures.
root=fileparts(fileparts(mfilename('fullpath')));
old=path;cleanup=onCleanup(@()path(old)); %#ok<NASGU>
addpath(fullfile(root,'YawAnalysis'));
calls={ ...
    @()yawSpatialLookAheadAnalysis(struct('saveFigures',true)), ...
    @()yawSpatialLookAheadMonteCarloAnalysis(struct('plotCombined',true)), ...
    @()yawSpatialLookAheadMonteCarloAnalysis(struct('individualDatasets',"lilyRun_20190309")), ...
    @()yawDistributionAnalysis(struct('saveFigures',true))};
for i=1:numel(calls)
    rejected=false;
    try
        calls{i}();
    catch err
        assert(strcmp(err.identifier,'release:nonManuscriptPlot'), ...
            'Unexpected error: %s',err.message);
        rejected=true;
    end
    assert(rejected,'Non-manuscript plotting option %d was accepted.',i);
end
fprintf('PASS four non-manuscript plot requests rejected before loading data.\n');
% Existing saved results contain all five pursuits. Empty selection must
% still draw only the two manuscript examples, including in plot-only mode.
source=fullfile(root,'outputs','yaw','spatial','MonteCarloResults', ...
    'yaw_spatial_lookahead_monte_carlo_analysis.mat');
assert(isfile(source),'Run run_saved_analysis before this check.');
out=tempname;mkdir(out);folderCleanup=onCleanup(@()rmdir(out,'s')); %#ok<NASGU>
copyfile(source,fullfile(out,'yaw_spatial_lookahead_monte_carlo_analysis.mat'));
yawSpatialLookAheadMonteCarloAnalysis(struct('plotOnly',true, ...
    'outputDir',out,'individualDatasets',strings(0,1),'figureResolution',72));
files=dir(fullfile(out,'*.png'));
expected=sort(["julesTopRun1_20171209_spatial_lookahead_monte_carlo_timeseries.png", ...
    "julesTopRun11_20170829_spatial_lookahead_monte_carlo_timeseries.png"]);
assert(isequal(sort(string({files.name})),expected),'Unexpected spatial figure exports.');
fprintf('PASS direct plot-only call exports exactly the two S15 examples.\n');
end
