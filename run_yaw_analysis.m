function run_yaw_analysis(sourceRoot, secondary)
%RUN_YAW_ANALYSIS Yaw geometry/states (18), ECDFs (7), signed errors (5).
% run_yaw_analysis([],true) also computes spatial lead and range closure.
if nargin<1, sourceRoot=[]; end
if nargin<2, secondary=false; end
C=prepare_analysis_inputs(sourceRoot,'yaw',true);
old=path; cleanup=onCleanup(@()path(old)); %#ok<NASGU>
addpath(fullfile(C.root,'YawAnalysis'));
out=fullfile(C.root,'outputs','yaw');
yawDistributionAnalysis(struct('configuration',C.configuration, ...
    'frameSelection',C.frames,'ecdfOnly',true,'saveFigures',false,'outputDir',out));
yawGeometryAndStates(C,fullfile(out,'geometry'));
if secondary
    dataRoot=fullfile(C.dataRoot,'yaw_tracking'); spatial=fullfile(out,'spatial');
    yawSpatialLookAheadAnalysis(struct('dataRoot',dataRoot,'outputDir',spatial,'extendEndpoint',true,'saveFigures',false));
    mcDir=fullfile(spatial,'MonteCarloResults');
    yawSpatialLookAheadMonteCarloAnalysis(struct('dataRoot',dataRoot, ...
        'nominalResultDir',spatial,'outputDir',mcDir,'extendEndpoint',true,'maxMC',500,'plotCombined',false, ...
        'individualDatasets',["julesTopRun1_20171209","julesTopRun11_20170829"], ...
        'individualTimeEndExclusive',struct('julesTopRun11_20170829',.64)));
    plotSpatialLeadRatioAcrossPursuits(struct('nominalCsv',fullfile(spatial,'yaw_spatial_lookahead_by_frame.csv'), ...
        'monteCarloMat',fullfile(mcDir,'yaw_spatial_lookahead_monte_carlo_analysis.mat'),'outputDir',spatial));
    handoff=fullfile(out,'head_torso');
    H=headTorsoHandoff(struct('dataRoot',dataRoot,'outputDir',handoff,'nMonteCarlo',500));
    plotHeadTorsoHandoff(H,struct('dataRoot',dataRoot,'resultDir',handoff));
end
fprintf('Yaw tables: %s. Render figures with python YawAnalysis/plot_yaw.py\n',out);
end
