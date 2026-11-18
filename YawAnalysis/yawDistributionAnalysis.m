function R = yawDistributionAnalysis(opts)
% yawDistributionAnalysis  Numerical source tables for yaw-error distributions.
%
% R = yawDistributionAnalysis(opts)
% Called by run_yaw_analysis with explicit cohort and frame IDs.
% Computes per-run and animal-equal ECDFs from paired MC-median yaw errors.

arguments
    opts struct = struct()
end

thisDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
opts = applyDefaults(opts, struct( ...
    'dataRoot', fullfile(repoRoot, 'outputs', 'analysis_inputs', 'yaw_tracking'), ...
    'outputDir', fullfile(repoRoot, 'outputs', 'yaw'), ...
    'configuration', table(), ...
    'frameSelection', table(), ...
    'ecdfOnly', false, ...
    'saveFigures', false, ...
    'referenceThresholdDeg', 10, ...
    'ecdfGridStepDeg', 0.5, ...
    'ecdfMaxDeg', 180, ...
    'signedGridStepDeg', 0.5, ...
    'signedXLimDeg', [-40 140], ...
    'signedBandwidthDeg', 5, ...
    'showFigures', false, ...
    'showTitle', false, ...
    'showSignedPanelTitles', true, ...
    'figureResolution', 600));

assert(~opts.saveFigures,'release:nonManuscriptPlot', ...
    'Render the manuscript cohort ECDFs with python YawAnalysis/plot_yaw.py.');

assert(opts.referenceThresholdDeg > 0, 'referenceThresholdDeg must be positive.');
assert(opts.ecdfGridStepDeg > 0, 'ecdfGridStepDeg must be positive.');
assert(opts.ecdfMaxDeg > 0, 'ecdfMaxDeg must be positive.');
assert(opts.signedGridStepDeg > 0, 'signedGridStepDeg must be positive.');
assert(opts.signedBandwidthDeg > 0, 'signedBandwidthDeg must be positive.');
assert(numel(opts.signedXLimDeg) == 2 && diff(opts.signedXLimDeg) > 0, ...
    'signedXLimDeg must be an increasing two-element vector.');

if ~exist(opts.outputDir, 'dir')
    mkdir(opts.outputDir);
end

if isempty(opts.configuration)
    error('Supply configuration and frameSelection via run_yaw_analysis.');
else
    C = opts.configuration;
    assert(istable(C) && all(ismember({'dataset','animal','context','dataDir'}, ...
        C.Properties.VariableNames)), 'configuration requires dataset, animal, context and dataDir.');
    assert(numel(unique(C.dataset))==height(C), 'Duplicate configured datasets.');
end
fprintf('\n==== yawDistributionAnalysis ====\n');
fprintf('Included high-speed pursuits: %d\n', height(C));
fprintf('Output directory: %s\n', opts.outputDir);

E = buildAbsoluteErrorEcdf(C, opts);
R = struct('configuration',C,'ecdf',E,'opts',opts);
save(fullfile(opts.outputDir,'yaw_absolute_error_analysis.mat'),'R');

fprintf('Saved yaw distribution analysis and figure-source tables.\n');
end

function E = buildAbsoluteErrorEcdf(C, opts)
gridDeg = (0:opts.ecdfGridStepDeg:opts.ecdfMaxDeg).';
nRuns = height(C);
headByRun = nan(numel(gridDeg), nRuns);
torsoByRun = nan(numel(gridDeg), nRuns);
runRows = cell(nRuns, 1);
summaryRows = cell(nRuns, 1);
frameRows = cell(nRuns, 1);

for i = 1:nRuns
    dataset = char(C.dataset(i));
    if ismember('dataDir', C.Properties.VariableNames)
        dataDir = char(C.dataDir(i));
    else
        dataDir = fullfile(opts.dataRoot, dataset);
    end
    mcPath = fullfile(dataDir, 'MonteCarloResults', ...
        ['mcMain_' dataset '.mat']);
    assert(isfile(mcPath), 'Missing Monte Carlo yaw result: %s', mcPath);
    M = load(mcPath, 'angleNames', 'medianDeg', 'absFrames');
    required = {'angleNames','medianDeg','absFrames'};
    for fi = 1:numel(required)
        assert(isfield(M, required{fi}), '%s lacks %s.', mcPath, required{fi});
    end

    gammaColumn = find(strcmp(M.angleNames, 'gamma_headYawErr'), 1);
    deltaColumn = find(strcmp(M.angleNames, 'delta_bodyYawErr'), 1);
    assert(~isempty(gammaColumn) && ~isempty(deltaColumn), ...
        '%s lacks gamma/delta angle columns.', mcPath);

    gammaMedian = double(M.medianDeg(:,gammaColumn));
    deltaMedian = double(M.medianDeg(:,deltaColumn));
    frames = double(M.absFrames(:));
    assert(numel(unique(frames))==numel(frames), 'Duplicate source frames: %s.', dataset);
    if ~isempty(opts.frameSelection)
        selected = double(opts.frameSelection.absFrame(opts.frameSelection.dataset==C.dataset(i)));
        assert(~isempty(selected) && numel(unique(selected))==numel(selected), ...
            'Empty or duplicate selected frames for %s.', dataset);
        [found, ix] = ismember(selected, frames);
        assert(all(found), 'Missing selected Monte Carlo frame for %s.', dataset);
        frames = frames(ix); gammaMedian = gammaMedian(ix); deltaMedian = deltaMedian(ix);
        assert(all(isfinite(gammaMedian)&isfinite(deltaMedian)), ...
            'Selected frames must have paired finite head and torso errors: %s.', dataset);
    end

    gammaStored = gammaMedian; deltaStored = deltaMedian;
    % Circular medians are stored as nominal + median deviation and can cross
    % the +/-180 boundary. Absolute yaw error is the shortest angular distance.
    gammaMedian = mod(gammaMedian+180,360)-180;
    deltaMedian = mod(deltaMedian+180,360)-180;
    headByRun(:,i) = empiricalCdf(abs(gammaMedian), gridDeg);
    torsoByRun(:,i) = empiricalCdf(abs(deltaMedian), gridDeg);

    runRows{i} = table(repmat(C.dataset(i), numel(gridDeg), 1), ...
        repmat(C.context(i), numel(gridDeg), 1), gridDeg, ...
        headByRun(:,i), torsoByRun(:,i), ...
        'VariableNames', {'dataset','context','absoluteYawDeg', ...
        'headPosteriorMedianTrajectoryCdf','torsoPosteriorMedianTrajectoryCdf'});

    summaryRows{i} = posteriorMedianSummaryRow(C.dataset(i), C.context(i), ...
        numel(frames), gammaMedian, deltaMedian, opts.referenceThresholdDeg);
    frameRows{i} = table(repmat(C.dataset(i), numel(frames), 1), ...
        repmat(C.animal(i), numel(frames), 1), repmat(C.context(i), numel(frames), 1), ...
        frames, gammaMedian, deltaMedian, abs(gammaMedian), abs(deltaMedian), gammaStored, deltaStored, ...
        'VariableNames', {'dataset','animal','context','absFrame','headMedianYawDeg', ...
        'torsoMedianYawDeg','headAbsoluteYawDeg','torsoAbsoluteYawDeg', ...
        'headStoredMedianYawDeg','torsoStoredMedianYawDeg'});
end

[headMean, headSE] = meanSeColumns(headByRun);
[torsoMean, torsoSE] = meanSeColumns(torsoByRun);

animals = unique(C.animal, 'stable');
headByAnimal = nan(numel(gridDeg), numel(animals));
torsoByAnimal = nan(numel(gridDeg), numel(animals));
for ai = 1:numel(animals)
    use = C.animal == animals(ai);
    headByAnimal(:,ai) = mean(headByRun(:,use), 2, 'omitnan');
    torsoByAnimal(:,ai) = mean(torsoByRun(:,use), 2, 'omitnan');
end
[headAnimalMean, headAnimalSE] = meanSeColumns(headByAnimal);
[torsoAnimalMean, torsoAnimalSE] = meanSeColumns(torsoByAnimal);

T = table(gridDeg, headMean, headSE, torsoMean, torsoSE, ...
    'VariableNames', {'absoluteYawDeg','headRunEqualCdfMean','headRunSE', ...
    'torsoRunEqualCdfMean','torsoRunSE'});
writetable(T, fullfile(opts.outputDir, 'yaw_absolute_error_ecdf.csv'));
writetable(table(gridDeg, headAnimalMean, headAnimalSE, torsoAnimalMean, torsoAnimalSE, ...
    'VariableNames', {'absoluteYawDeg','headAnimalEqualCdfMean','headAnimalSE', ...
    'torsoAnimalEqualCdfMean','torsoAnimalSE'}), ...
    fullfile(opts.outputDir, 'yaw_absolute_error_ecdf_animal_equal.csv'));
writetable(vertcat(runRows{:}), ...
    fullfile(opts.outputDir, 'yaw_absolute_error_ecdf_by_run.csv'));
writetable(vertcat(summaryRows{:}), ...
    fullfile(opts.outputDir, 'yaw_absolute_error_summary_by_run.csv'));
writetable(vertcat(frameRows{:}), ...
    fullfile(opts.outputDir, 'yaw_absolute_error_by_frame.csv'));
writetable(C, fullfile(opts.outputDir, 'yaw_absolute_error_cohort.csv'));

E = struct('table', T, 'gridDeg', gridDeg, 'headByRun', headByRun, ...
    'torsoByRun', torsoByRun, 'headMean', headMean, 'headSE', headSE, ...
    'torsoMean', torsoMean, 'torsoSE', torsoSE, 'animals', animals, ...
    'headByAnimal', headByAnimal, 'torsoByAnimal', torsoByAnimal, ...
    'headAnimalMean', headAnimalMean, 'headAnimalSE', headAnimalSE, ...
    'torsoAnimalMean', torsoAnimalMean, 'torsoAnimalSE', torsoAnimalSE, ...
    'frameTable', vertcat(frameRows{:}), 'summaryByRun', vertcat(summaryRows{:}));
end

function F = empiricalCdf(values, grid)
values = values(isfinite(values));
if isempty(values)
    F = nan(size(grid));
else
    F = arrayfun(@(g) mean(values <= g), grid);
end
F = F(:);
end

function [mu, se] = meanSeColumns(X)
mu = mean(X, 2, 'omitnan');
n = sum(isfinite(X), 2);
se = std(X, 0, 2, 'omitnan') ./ sqrt(n);
se(n < 2) = NaN;
end

function row = posteriorMedianSummaryRow(dataset, context, nFrames, gamma, delta, threshold)
gamma = abs(double(gamma(:)));
delta = abs(double(delta(:)));
row = table(dataset, context, nFrames, ...
    median(gamma, 'omitnan'), median(delta, 'omitnan'), ...
    percentile(gamma, 90), percentile(delta, 90), ...
    mean(gamma <= threshold, 'omitnan'), mean(delta <= threshold, 'omitnan'), ...
    'VariableNames', {'dataset','context','nFrames','headMedianAbsDeg', ...
    'torsoMedianAbsDeg','headQ90AbsDeg','torsoQ90AbsDeg', ...
    'headWithinReference','torsoWithinReference'});
end

function q = percentile(x, p)
x = sort(x(isfinite(x)));
if isempty(x)
    q = NaN;
    return;
end
if isscalar(x)
    q = x;
    return;
end
position = 1 + (numel(x)-1) * p/100;
lo = floor(position); hi = ceil(position);
q = x(lo) + (position-lo) * (x(hi)-x(lo));
end

function opts = applyDefaults(opts, defaults)
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end

