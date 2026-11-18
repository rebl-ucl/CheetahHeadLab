function P = pitchTightnessAnalysis(datasetList, opts)
% pitchTightnessAnalysis  Pitch-error ECDF and run-level tightness.
%
% P = pitchTightnessAnalysis(datasetList, opts)
%
% Inputs
%   datasetList  Optional cell/string list of "group/dataset" paths beneath
%                dataRoot. If omitted, the shared pitch-analysis cohort is used.
%   opts         Optional struct; see runPitchAnalysis.m.
%
% Outputs (CSV, PNG, PDF, and pitch_tightness_analysis.mat)
%   pitch_absolute_error_ecdf.*
%       Individual-sequence ECDFs plus an animal- or run-equal mean ECDF of
%       |alpha| and |beta|. The optional 10-degree line is a visual reference,
%       not part of the ECDF calculation.
%   pitch_centered_error_distribution.*
%       Run-centred signed-error distributions with the same run- or
%       animal-equal aggregation convention.
% alpha = theta_h_pitch - lambda_h_pitch (head-to-lure pitch error)
% beta  = theta_t_pitch - lambda_t_pitch (torso-to-lure pitch error)

arguments
    datasetList = []
    opts struct = struct()
end

% Make the interactive no-argument call useful while preserving the original
% behaviour when a dataset list is supplied explicitly.
if isempty(datasetList)
    C = pitchAnalysisConfig();
    datasetList = C.datasets;
    if isempty(fieldnames(opts))
        opts = C.commonOpts;
    end
end

thisDir = fileparts(mfilename('fullpath'));
lureDir = thisDir;
repoRoot = fileparts(thisDir);
opts = applyDefaults(opts, struct( ...
    'dataRoot', fullfile(repoRoot, 'outputs', 'analysis_inputs'), ...
    'outputDir', fullfile(repoRoot, 'outputs', 'pitch', 'tightness_run_equal'), ...
    'resultFileName', 'lureAnalysis_results_trajectory.mat', ...
    'frameRanges', {cell(0, 2)}, ...
    'frameRangeMode', 'absFrame', ...
    'excludeInterpolated', false, ...
    'useDespiked', false, ...
    'referenceThresholdDeg', 10, ...
    'ecdfGridStepDeg', 0.25, ...
    'ecdfMaxDeg', [], ...
    'ecdfWeighting', 'animal', ...
    'saveCenteredDistribution', true, ...
    'centeredBinWidthDeg', 5, ...
    'centeredMaxDeg', [], ...
    'errorBarPositionsDeg', -15:5:15, ...
    'counterRotationLagSec', 0.05, ...
    'counterRotationSmoothSec', 0.05, ...
    'minTorsoChangeDeg', 1.0, ...
    'minCounterRotationPairs', 10, ...
    'maxScatterPointsPerRun', 180, ...
    'saveCounterRotationOutputs', false, ...
    'saveFigures', true, ...
    'showFigures', false, ...
    'showFigureTitles', false, ...
    'titleFontSize', 16, ...
    'figureResolution', 600));

opts.dataRoot = char(opts.dataRoot);
opts.outputDir = char(opts.outputDir);
opts.resultFileName = char(opts.resultFileName);
opts.frameRangeMode = validatestring(opts.frameRangeMode, {'absFrame','alignedRow'});
opts.ecdfWeighting = validatestring(opts.ecdfWeighting, {'animal','run'});
assert(opts.ecdfGridStepDeg > 0, 'ecdfGridStepDeg must be positive.');
assert(opts.centeredBinWidthDeg > 0, 'centeredBinWidthDeg must be positive.');
assert(opts.referenceThresholdDeg > 0, 'referenceThresholdDeg must be positive.');
assert(opts.counterRotationLagSec > 0, 'counterRotationLagSec must be positive.');
assert(opts.counterRotationSmoothSec >= 0, 'counterRotationSmoothSec cannot be negative.');
assert(opts.minTorsoChangeDeg >= 0, 'minTorsoChangeDeg cannot be negative.');

if ~exist(opts.outputDir, 'dir')
    mkdir(opts.outputDir);
end

datasetList = normaliseList(datasetList);
assert(~isempty(datasetList), 'datasetList cannot be empty.');

runData = struct('key', {}, 'group', {}, 'dataset', {}, 'animal', {}, ...
    'alphaDeg', {}, 'betaDeg', {}, 'absFrames', {}, ...
    'requiredRelDeg', {}, 'actualRelDeg', {}, ...
    'dTorsoDeg', {}, 'dNeckDeg', {}, 'dHeadDeg', {});
summaryRows = {};
frameRows = {};
counterRows = {};

fprintf('\n==== pitchTightnessAnalysis ====' );
fprintf('\nResult file: %s\n', opts.resultFileName);

for i = 1:numel(datasetList)
    key = normaliseKey(datasetList{i});
    [group, dataset] = splitDatasetKey(key);
    resultFile = fullfile(opts.dataRoot, group, dataset, opts.resultFileName);
    if ~isfile(resultFile)
        error('pitchTightnessAnalysis:missingResult', ...
            'Missing %s for %s.', opts.resultFileName, key);
    end

    S = load(resultFile, 'runResults');
    assert(isfield(S, 'runResults'), 'runResults is absent from %s.', resultFile);
    R = S.runResults;
    requiredFields = {'alpha','beta','theta_h_pitch','theta_t_pitch', ...
        'lambda_h_pitch','pitchAnalysisMask'};
    for fi = 1:numel(requiredFields)
        assert(isfield(R, requiredFields{fi}), ...
            '%s lacks runResults.%s. Re-run lureAnalysis.', key, requiredFields{fi});
    end

    if opts.useDespiked && isfield(R, 'alphaF') && isfield(R, 'betaF')
        alpha = R.alphaF(:);
        beta = R.betaF(:);
        seriesSource = "despiked";
    else
        if opts.useDespiked
            warning('pitchTightnessAnalysis:noDespiked', ...
                '%s has no alphaF/betaF; using raw errors.', key);
        end
        alpha = R.alpha(:);
        beta = R.beta(:);
        seriesSource = "raw";
    end

    thetaH = R.theta_h_pitch(:);
    thetaT = R.theta_t_pitch(:);
    lambdaH = R.lambda_h_pitch(:);
    mask = logical(R.pitchAnalysisMask(:));
    n = min([numel(alpha), numel(beta), numel(thetaH), numel(thetaT), ...
        numel(lambdaH), numel(mask)]);
    alpha = alpha(1:n); beta = beta(1:n);
    thetaH = thetaH(1:n); thetaT = thetaT(1:n); lambdaH = lambdaH(1:n);
    mask = mask(1:n);
    absFrames = absoluteFrames(R, n);

    if opts.excludeInterpolated && isfield(R, 'isInterp')
        interp = false(n, 1);
        tmp = logical(R.isInterp(:));
        interp(1:min(n, numel(tmp))) = tmp(1:min(n, numel(tmp)));
        mask = mask & ~interp;
    end
    ranges = lookupRanges(key, dataset, opts.frameRanges);
    if ~isempty(ranges)
        if strcmp(opts.frameRangeMode, 'absFrame')
            rangeValues = absFrames;
        else
            rangeValues = (1:n).';
        end
        mask = mask & inRanges(rangeValues, ranges);
    end
    mask = mask & isfinite(alpha) & isfinite(beta) & isfinite(thetaH) & ...
        isfinite(thetaT) & isfinite(lambdaH);

    alphaDegAll = rad2deg(alpha);
    betaDegAll = rad2deg(beta);
    alphaDeg = alphaDegAll(mask);
    betaDeg = betaDegAll(mask);
    keptFrames = absFrames(mask);
    if isempty(alphaDeg)
        warning('pitchTightnessAnalysis:noFrames', '%s contributed no frames.', key);
        continue;
    end

    animal = inferAnimal(dataset);
    nFrames = numel(alphaDeg);
    headMedianAbsDeg = median(abs(alphaDeg));
    bodyMedianAbsDeg = median(abs(betaDeg));
    headQ80Deg = percentile(abs(alphaDeg), 80);
    bodyQ80Deg = percentile(abs(betaDeg), 80);
    headQ90Deg = percentile(abs(alphaDeg), 90);
    bodyQ90Deg = percentile(abs(betaDeg), 90);
    headRmsDeg = sqrt(mean(alphaDeg.^2));
    bodyRmsDeg = sqrt(mean(betaDeg.^2));
    headWithinReference = mean(abs(alphaDeg) <= opts.referenceThresholdDeg);
    bodyWithinReference = mean(abs(betaDeg) <= opts.referenceThresholdDeg);
    headSignedMedianDeg = median(alphaDeg);
    bodySignedMedianDeg = median(betaDeg);

    % Pitch-plane counter-rotation. Each consecutive valid segment is
    % smoothed separately, so neither excluded gaps nor curated ranges are
    % bridged. actualRel-requiredRel equals head-lure pitch error before the
    % optional smoothing, by construction.
    fs = getScalarField(R, 'fs', 1);
    lagFrames = max(1, round(opts.counterRotationLagSec * fs));
    smoothFrames = smoothingWindow(opts.counterRotationSmoothSec, fs);
    segments = consecutiveSegments(find(mask), absFrames);
    requiredRelDeg = [];
    actualRelDeg = [];
    dTorsoDeg = [];
    dNeckDeg = [];
    dHeadDeg = [];
    counterStartFrame = [];
    counterEndFrame = [];

    for si = 1:numel(segments)
        idx = segments{si};
        if isempty(idx), continue; end
        h = smoothAngle(thetaH(idx), smoothFrames);
        t = smoothAngle(thetaT(idx), smoothFrames);
        l = smoothAngle(lambdaH(idx), smoothFrames);
        actual = rad2deg(h - t);
        required = rad2deg(l - t);
        actualRelDeg = [actualRelDeg; actual(:)]; %#ok<AGROW>
        requiredRelDeg = [requiredRelDeg; required(:)]; %#ok<AGROW>

        if numel(idx) <= lagFrames, continue; end
        a = 1:(numel(idx) - lagFrames);
        b = a + lagFrames;
        dT = rad2deg(t(b) - t(a));
        dN = rad2deg((h(b) - t(b)) - (h(a) - t(a)));
        dH = rad2deg(h(b) - h(a));
        good = isfinite(dT) & isfinite(dN) & isfinite(dH) & ...
            abs(dT) >= opts.minTorsoChangeDeg;
        dTorsoDeg = [dTorsoDeg; dT(good(:))]; %#ok<AGROW>
        dNeckDeg = [dNeckDeg; dN(good(:))]; %#ok<AGROW>
        dHeadDeg = [dHeadDeg; dH(good(:))]; %#ok<AGROW>
        counterStartFrame = [counterStartFrame; absFrames(idx(a(good)))]; %#ok<AGROW>
        counterEndFrame = [counterEndFrame; absFrames(idx(b(good)))]; %#ok<AGROW>
    end

    [counterSlope, counterGain, counterCorrelation, opposingFraction] = ...
        counterMetrics(dTorsoDeg, dNeckDeg, opts.minCounterRotationPairs);
    [requiredSlope, requiredR2, requiredRmseDeg] = ...
        requiredCompensationMetrics(requiredRelDeg, actualRelDeg);

    summaryRows{end+1} = table(string(key), string(group), string(dataset), ... %#ok<AGROW>
        string(animal), string(seriesSource), string(rangesText(ranges)), ...
        nFrames, fs, headSignedMedianDeg, bodySignedMedianDeg, ...
        headMedianAbsDeg, bodyMedianAbsDeg, headQ80Deg, bodyQ80Deg, ...
        headQ90Deg, bodyQ90Deg, headRmsDeg, bodyRmsDeg, ...
        headWithinReference, bodyWithinReference, numel(dTorsoDeg), ...
        lagFrames / fs, smoothFrames / fs, counterSlope, counterGain, ...
        counterCorrelation, opposingFraction, requiredSlope, requiredR2, requiredRmseDeg, ...
        'VariableNames', {'datasetKey','group','dataset','animal','seriesSource', ...
        'frameRanges','nFrames','fs','headSignedMedianDeg','torsoSignedMedianDeg', ...
        'headMedianAbsDeg','torsoMedianAbsDeg','headQ80AbsDeg','torsoQ80AbsDeg', ...
        'headQ90AbsDeg','torsoQ90AbsDeg','headRmsDeg','torsoRmsDeg', ...
        'headWithinReference','torsoWithinReference','nCounterPairs', ...
        'counterLagSecActual','counterSmoothSecActual','counterSlope', ...
        'counterRotationGain','counterCorrelation','opposingFraction', ...
        'requiredCompensationSlope','requiredCompensationR2','requiredCompensationRmseDeg'});

    frameRows{end+1} = table(repmat(string(key), nFrames, 1), ... %#ok<AGROW>
        repmat(string(animal), nFrames, 1), keptFrames(:), alphaDeg(:), betaDeg(:), ...
        abs(alphaDeg(:)), abs(betaDeg(:)), ...
        'VariableNames', {'datasetKey','animal','absFrame','alphaDeg','betaDeg', ...
        'absAlphaDeg','absBetaDeg'});

    if ~isempty(dTorsoDeg)
        counterRows{end+1} = table(repmat(string(key), numel(dTorsoDeg), 1), ... %#ok<AGROW>
            repmat(string(animal), numel(dTorsoDeg), 1), counterStartFrame(:), ...
            counterEndFrame(:), dTorsoDeg(:), dNeckDeg(:), dHeadDeg(:), ...
            'VariableNames', {'datasetKey','animal','startAbsFrame','endAbsFrame', ...
            'deltaTorsoPitchDeg','deltaHeadRelativeTorsoPitchDeg','deltaHeadPitchDeg'});
    end

    runData(end+1) = struct('key', string(key), 'group', string(group), ... %#ok<AGROW>
        'dataset', string(dataset), 'animal', string(animal), ...
        'alphaDeg', alphaDeg(:), 'betaDeg', betaDeg(:), 'absFrames', keptFrames(:), ...
        'requiredRelDeg', requiredRelDeg(:), 'actualRelDeg', actualRelDeg(:), ...
        'dTorsoDeg', dTorsoDeg(:), 'dNeckDeg', dNeckDeg(:), 'dHeadDeg', dHeadDeg(:));

    fprintf('  %-52s n=%4d | q90 head %6.2f deg, torso %6.2f deg\n', ...
        key, nFrames, headQ90Deg, bodyQ90Deg);
end

runSummary = vertcatSafe(summaryRows);
frameTable = vertcatSafe(frameRows);
counterTable = vertcatSafe(counterRows);
assert(~isempty(runSummary), 'None of the requested datasets contributed pitch frames.');

if ~opts.saveCounterRotationOutputs
    counterTable = table();
    exploratoryColumns = {'nCounterPairs','counterLagSecActual', ...
        'counterSmoothSecActual','counterSlope','counterRotationGain', ...
        'counterCorrelation','opposingFraction','requiredCompensationSlope', ...
        'requiredCompensationR2','requiredCompensationRmseDeg'};
    runSummary(:, exploratoryColumns) = [];
end

writetable(runSummary, fullfile(opts.outputDir, 'pitch_run_level_tightness.csv'));
writetable(frameTable, fullfile(opts.outputDir, 'pitch_frame_level_errors.csv'));
if opts.saveCounterRotationOutputs && ~isempty(counterTable)
    writetable(counterTable, fullfile(opts.outputDir, 'pitch_head_neck_counterrotation_pairs.csv'));
end

ecdf = plotAbsoluteErrorEcdf(runData, opts);
if opts.saveCenteredDistribution
    centered = plotCenteredErrorDistribution(runData, opts);
else
    centered = struct([]);
end
% Optional counter-rotation tables are retained without a legacy figure.

P = struct();
P.opts = opts;
P.runSummary = runSummary;
P.frameTable = frameTable;
P.counterTable = counterTable;
P.ecdf = ecdf;
P.centered = centered;
P.outputDir = opts.outputDir;
save(fullfile(opts.outputDir, 'pitch_tightness_analysis.mat'), 'P');
writeInterpretationNotes(P);
end

function E = plotAbsoluteErrorEcdf(runData, opts)
pooledBody = vertcat(runData.betaDeg);
pooledHead = vertcat(runData.alphaDeg);
if isempty(opts.ecdfMaxDeg)
    xMax = max(30, percentile(abs([pooledHead; pooledBody]), 95));
    xMax = ceil(xMax / 5) * 5;
else
    xMax = opts.ecdfMaxDeg;
end
gridDeg = (0:opts.ecdfGridStepDeg:xMax).';
nRuns = numel(runData);
headByRun = nan(numel(gridDeg), nRuns);
bodyByRun = nan(numel(gridDeg), nRuns);
runKeys = strings(nRuns, 1);
animals = strings(nRuns, 1);
for i = 1:nRuns
    headByRun(:,i) = empiricalCdf(abs(runData(i).alphaDeg), gridDeg);
    bodyByRun(:,i) = empiricalCdf(abs(runData(i).betaDeg), gridDeg);
    runKeys(i) = runData(i).key;
    animals(i) = runData(i).animal;
end

if strcmp(opts.ecdfWeighting, 'animal')
    units = unique(animals, 'stable');
    headByUnit = nan(numel(gridDeg), numel(units));
    bodyByUnit = nan(numel(gridDeg), numel(units));
    for i = 1:numel(units)
        m = animals == units(i);
        headByUnit(:,i) = mean(headByRun(:,m), 2, 'omitnan');
        bodyByUnit(:,i) = mean(bodyByRun(:,m), 2, 'omitnan');
    end
else
    units = runKeys;
    headByUnit = headByRun;
    bodyByUnit = bodyByRun;
end
[headMean, headSE] = meanSeColumns(headByUnit);
[bodyMean, bodySE] = meanSeColumns(bodyByUnit);

T = table(gridDeg, headMean, headSE, bodyMean, bodySE, ...
    'VariableNames', {'absoluteErrorDeg','headMeanCdf','headSE','torsoMeanCdf','torsoSE'});
writetable(T, fullfile(opts.outputDir, 'pitch_absolute_error_ecdf.csv'));

longRows = {};
for i = 1:nRuns
    longRows{end+1} = table(repmat(runKeys(i), numel(gridDeg), 1), ... %#ok<AGROW>
        repmat(animals(i), numel(gridDeg), 1), gridDeg, headByRun(:,i), bodyByRun(:,i), ...
        'VariableNames', {'datasetKey','animal','absoluteErrorDeg','headCdf','torsoCdf'});
end
writetable(vertcatSafe(longRows), fullfile(opts.outputDir, 'pitch_absolute_error_ecdf_by_run.csv'));

[headColor, torsoColor] = pitchColors();
visible = onoff(opts.showFigures);
f = figure('Color', 'w', 'Visible', visible, 'Position', [100 100 760 560]);
hold on;
for i = 1:nRuns
    plot(gridDeg, bodyByRun(:,i), '-', 'Color', mixWithWhite(torsoColor, 0.78), ...
        'LineWidth', 0.7, 'HandleVisibility', 'off');
    plot(gridDeg, headByRun(:,i), '-', 'Color', mixWithWhite(headColor, 0.78), ...
        'LineWidth', 0.7, 'HandleVisibility', 'off');
end
shadeBand(gridDeg, bodyMean-bodySE, bodyMean+bodySE, torsoColor);
shadeBand(gridDeg, headMean-headSE, headMean+headSE, headColor);
pBody = plot(gridDeg, bodyMean, '-', 'Color', torsoColor, 'LineWidth', 2.5, ...
    'DisplayName', 'Torso-lure pitch error');
pHead = plot(gridDeg, headMean, '-', 'Color', headColor, 'LineWidth', 2.5, ...
    'DisplayName', 'Head-lure pitch error');
pBodyBand = patch(nan, nan, torsoColor, 'FaceAlpha', 0.13, ...
    'EdgeColor', 'none', 'DisplayName', 'Torso: +/- 1 SE');
pHeadBand = patch(nan, nan, headColor, 'FaceAlpha', 0.13, ...
    'EdgeColor', 'none', 'DisplayName', 'Head: +/- 1 SE');
xline(opts.referenceThresholdDeg, '--k', sprintf('%g deg', opts.referenceThresholdDeg), ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', ...
    'HandleVisibility', 'off');
hold off; box off; grid on;
xlim([0 xMax]); ylim([0 1]);
xlabel('Absolute pitch error (deg)');
ylabel('Cumulative proportion of frames');
legend([pHead pHeadBand pBody pBodyBand], 'Location', 'southeast', 'Box', 'on');
if opts.showFigureTitles
    title(sprintf('Absolute pitch error distributions (%s-equal mean)', ...
        opts.ecdfWeighting), 'FontSize', opts.titleFontSize, 'FontWeight', 'bold');
end
saveFigurePair(f, fullfile(opts.outputDir, 'pitch_absolute_error_ecdf'), opts);

E = struct('table', T, 'gridDeg', gridDeg, 'headByRun', headByRun, ...
    'bodyByRun', bodyByRun, 'headByUnit', headByUnit, 'bodyByUnit', bodyByUnit, ...
    'units', units, 'weighting', string(opts.ecdfWeighting));
end

function C = plotCenteredErrorDistribution(runData, opts)
% Run-centered signed pitch-error distributions (head vs torso), styled after
% Eitan et al. 2019 Proc. R. Soc. B figure 1c. Each run's signed error is
% shifted to its own mean before binning, so the pooled curve reflects
% within-run jitter (precision) and NOT the between-run pointing offset. This
% is a precision/stabilisation figure, not an accuracy figure.
nRuns = numel(runData);
headCentered = cell(nRuns, 1);
bodyCentered = cell(nRuns, 1);
for i = 1:nRuns
    a = runData(i).alphaDeg(:); a = a(isfinite(a));
    b = runData(i).betaDeg(:);  b = b(isfinite(b));
    headCentered{i} = a - mean(a);
    bodyCentered{i} = b - mean(b);
end

pooledAbs = abs([vertcat(headCentered{:}); vertcat(bodyCentered{:})]);
if isempty(opts.centeredMaxDeg)
    xMax = max(15, percentile(pooledAbs, 99));
    xMax = ceil(xMax / opts.centeredBinWidthDeg) * opts.centeredBinWidthDeg;
else
    xMax = opts.centeredMaxDeg;
end

% Bin edges placed so a bin is centred on 0.
bw = opts.centeredBinWidthDeg;
edges = (-xMax - bw/2) : bw : (xMax + bw/2);
centers = edges(1:end-1) + bw/2;

headByRun = nan(numel(centers), nRuns);
bodyByRun = nan(numel(centers), nRuns);
runKeys = strings(nRuns, 1);
animals = strings(nRuns, 1);
for i = 1:nRuns
    headByRun(:,i) = binProportion(headCentered{i}, edges);
    bodyByRun(:,i) = binProportion(bodyCentered{i}, edges);
    runKeys(i) = runData(i).key;
    animals(i) = runData(i).animal;
end

% Match the ECDF weighting convention: in animal-equal mode, first average
% run-level distributions within each animal, then average across animals.
if strcmp(opts.ecdfWeighting, 'animal')
    units = unique(animals, 'stable');
    headByUnit = nan(numel(centers), numel(units));
    bodyByUnit = nan(numel(centers), numel(units));
    for i = 1:numel(units)
        m = animals == units(i);
        headByUnit(:,i) = mean(headByRun(:,m), 2, 'omitnan');
        bodyByUnit(:,i) = mean(bodyByRun(:,m), 2, 'omitnan');
    end
else
    units = runKeys;
    headByUnit = headByRun;
    bodyByUnit = bodyByRun;
end

[headMean, headSE] = meanSeColumns(headByUnit);
[bodyMean, bodySE] = meanSeColumns(bodyByUnit);
[headL, headR, headFwhm] = halfMaxWidth(centers, headMean);
[bodyL, bodyR, bodyFwhm] = halfMaxWidth(centers, bodyMean);

T = table(centers(:), headMean, headSE, bodyMean, bodySE, ...
    'VariableNames', {'centeredErrorDeg','headMeanProp','headSE','torsoMeanProp','torsoSE'});
writetable(T, fullfile(opts.outputDir, 'pitch_centered_error_distribution.csv'));

[headColor, torsoColor] = pitchColors();
visible = onoff(opts.showFigures);
f = figure('Color', 'w', 'Visible', visible, 'Position', [100 100 760 560]);
hold on;
for i = 1:nRuns
    plot(centers, bodyByRun(:,i), '-', 'Color', mixWithWhite(torsoColor, 0.80), ...
        'LineWidth', 0.6, 'HandleVisibility', 'off');
    plot(centers, headByRun(:,i), '-', 'Color', mixWithWhite(headColor, 0.80), ...
        'LineWidth', 0.6, 'HandleVisibility', 'off');
end
shadeProportionBand(centers, bodyMean-bodySE, bodyMean+bodySE, torsoColor);
shadeProportionBand(centers, headMean-headSE, headMean+headSE, headColor);
pBody = plot(centers, bodyMean, '-', 'Color', torsoColor, 'LineWidth', 2.5, ...
    'DisplayName', 'Torso-lure pitch error');
pHead = plot(centers, headMean, '-', 'Color', headColor, 'LineWidth', 2.5, ...
    'DisplayName', 'Head-lure pitch error');
xline(0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
% Width-at-half-height markers.
if isfinite(headFwhm)
    plot([headL headR], (max(headMean)/2)*[1 1], '-', 'Color', headColor, ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
end
if isfinite(bodyFwhm)
    plot([bodyL bodyR], (max(bodyMean)/2)*[1 1], '-', 'Color', torsoColor, ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
end
% Discrete SE error bars at "stress" positions, matching Eitan et al. 1c.
ebPos = opts.errorBarPositionsDeg(:).';
for p = ebPos
    [dmin, bi] = min(abs(centers(:) - p));
    if dmin <= bw/2 + 1e-6
        errorbar(centers(bi), headMean(bi), headSE(bi), 'o', 'Color', headColor, ...
            'MarkerFaceColor', headColor, 'MarkerEdgeColor', 'w', 'MarkerSize', 6, ...
            'LineWidth', 1.4, 'CapSize', 6, 'HandleVisibility', 'off');
        errorbar(centers(bi), bodyMean(bi), bodySE(bi), 'o', 'Color', torsoColor, ...
            'MarkerFaceColor', torsoColor, 'MarkerEdgeColor', 'w', 'MarkerSize', 6, ...
            'LineWidth', 1.4, 'CapSize', 6, 'HandleVisibility', 'off');
    end
end
hold off; box off; grid on;
set(gca, 'XTick', -xMax:5:xMax);
xlim([-(xMax+bw) (xMax+bw)]);
xlabel('Run-centered pitch error (deg)');
ylabel('Proportion of frames');
legend([pHead pBody], 'Location', 'northeast');
if opts.showFigureTitles
    title(sprintf('Run-centered pitch-error distributions (%s-equal mean)', ...
        opts.ecdfWeighting), 'FontSize', opts.titleFontSize, ...
        'FontWeight', 'bold');
end
saveFigurePair(f, fullfile(opts.outputDir, 'pitch_centered_error_distribution'), opts);

C = struct('table', T, 'centers', centers(:), 'headByRun', headByRun, ...
    'bodyByRun', bodyByRun, 'headByUnit', headByUnit, ...
    'bodyByUnit', bodyByUnit, 'runKeys', runKeys, 'units', units, ...
    'weighting', string(opts.ecdfWeighting), 'binWidthDeg', bw, ...
    'headFwhmDeg', headFwhm, 'torsoFwhmDeg', bodyFwhm);
end

function p = binProportion(x, edges)
x = x(isfinite(x));
counts = histcounts(x, edges);
total = sum(counts);
if total == 0
    p = nan(numel(edges)-1, 1);
else
    p = counts(:) / total;
end
end

function [xL, xR, w] = halfMaxWidth(x, y)
x = x(:); y = y(:);
xL = NaN; xR = NaN; w = NaN;
valid = isfinite(y);
if nnz(valid) < 2, return; end
[ymax, pk] = max(y);
if ~isfinite(ymax) || ymax <= 0, return; end
half = ymax / 2;
xL = x(1);
for i = pk:-1:2
    if y(i) >= half && y(i-1) < half
        xL = interp1([y(i-1) y(i)], [x(i-1) x(i)], half);
        break;
    end
end
xR = x(end);
for i = pk:numel(y)-1
    if y(i) >= half && y(i+1) < half
        xR = interp1([y(i) y(i+1)], [x(i) x(i+1)], half);
        break;
    end
end
w = xR - xL;
end

function shadeProportionBand(x, lo, hi, color)
x = x(:); lo = max(0, lo(:)); hi = hi(:);
valid = isfinite(x) & isfinite(lo) & isfinite(hi);
if ~any(valid), return; end
fill([x(valid); flipud(x(valid))], [lo(valid); flipud(hi(valid))], color, ...
    'FaceAlpha', 0.13, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function [slope, gain, rho, opposing] = counterMetrics(dTorso, dNeck, minPairs)
slope = NaN; gain = NaN; rho = NaN; opposing = NaN;
valid = isfinite(dTorso) & isfinite(dNeck);
x = dTorso(valid); y = dNeck(valid);
if numel(x) < minPairs || std(x) == 0, return; end
b = [ones(numel(x),1), x(:)] \ y(:);
slope = b(2);
gain = -slope;
C = corrcoef(x, y);
if all(size(C) == [2 2]), rho = C(1,2); end
opposing = mean(x .* y < 0);
end

function [slope, r2, rmseDeg] = requiredCompensationMetrics(required, actual)
slope = NaN; r2 = NaN; rmseDeg = NaN;
valid = isfinite(required) & isfinite(actual);
x = required(valid); y = actual(valid);
if numel(x) < 3 || std(x) == 0, return; end
b = [ones(numel(x),1), x(:)] \ y(:);
yhat = [ones(numel(x),1), x(:)] * b;
slope = b(2);
ssRes = sum((y(:)-yhat).^2);
ssTot = sum((y(:)-mean(y)).^2);
if ssTot > 0, r2 = 1 - ssRes/ssTot; end
rmseDeg = sqrt(mean((y(:)-x(:)).^2));
end

function segments = consecutiveSegments(indices, absFrames)
segments = {};
if isempty(indices), return; end
cut = [true; diff(indices(:)) ~= 1 | diff(absFrames(indices)) ~= 1];
starts = find(cut);
ends = [starts(2:end)-1; numel(indices)];
for i = 1:numel(starts)
    segments{end+1} = indices(starts(i):ends(i)); %#ok<AGROW>
end
end

function x = smoothAngle(x, window)
x = unwrap(x(:));
if window > 1 && numel(x) >= 3
    w = min(window, numel(x));
    if mod(w,2)==0, w = w-1; end
    if w >= 3, x = movmedian(x, w); end
end
end

function w = smoothingWindow(sec, fs)
if sec <= 0
    w = 1;
    return;
end
w = max(3, round(sec * fs));
if mod(w,2)==0, w = w+1; end
end

function y = empiricalCdf(x, grid)
x = x(isfinite(x));
y = nan(size(grid));
if isempty(x), return; end
for i = 1:numel(grid)
    y(i) = mean(x <= grid(i));
end
end

function [mu, se] = meanSeColumns(X)
mu = mean(X, 2, 'omitnan');
n = sum(isfinite(X), 2);
se = std(X, 0, 2, 'omitnan') ./ sqrt(n);
se(n < 2) = NaN;
end

function shadeBand(x, lo, hi, color)
lo = max(0, lo); hi = min(1, hi);
valid = isfinite(x) & isfinite(lo) & isfinite(hi);
if ~any(valid), return; end
fill([x(valid); flipud(x(valid))], [lo(valid); flipud(hi(valid))], color, ...
    'FaceAlpha', 0.13, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function [headColor, torsoColor] = pitchColors()
headColor = [0.38 0.40 0.88];
torsoColor = [0.95 0.52 0.28];
end

function c = mixWithWhite(c, amount)
c = c*(1-amount) + amount;
end

function lim = symmetricLimit(x, y)
v = abs([x(:); y(:)]);
v = v(isfinite(v));
if isempty(v), m = 1; else, m = percentile(v, 98); end
m = max(m, 1);
lim = [-m m];
end

function idx = evenSample(n, maxN)
if n == 0
    idx = [];
elseif n <= maxN
    idx = (1:n).';
else
    idx = unique(round(linspace(1, n, maxN))).';
end
end

function saveFigurePair(f, basePath, opts)
if opts.saveFigures
    exportgraphics(f, [basePath '.png'], 'Resolution', opts.figureResolution);
    exportgraphics(f, [basePath '.pdf'], 'ContentType', 'vector');
end
if ~opts.showFigures, close(f); end
end

function writeInterpretationNotes(P)
path = fullfile(P.outputDir, 'PITCH_TIGHTNESS_NOTES.txt');
fid = fopen(path, 'w');
if fid < 0, return; end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, 'Pitch tightness analysis\n');
fprintf(fid, '========================\n\n');
fprintf(fid, 'Sequences: %d; animals: %d.\n', height(P.runSummary), ...
    numel(unique(P.runSummary.animal)));
fprintf(fid, 'ECDF weighting: %s-equal.\n', P.opts.ecdfWeighting);
fprintf(fid, 'Reference line: +/- %.1f deg (display only; ECDF and quantiles are threshold-free).\n\n', ...
    P.opts.referenceThresholdDeg);
fprintf(fid, ['Retained outputs quantify target-relative pitch tightness only. ' ...
    'Mechanistic counter-rotation is tested separately with exact 3-D and torso-defined event analyses.\n']);
end

function ranges = lookupRanges(key, dataset, entries)
ranges = [];
if isempty(entries), return; end
assert(iscell(entries) && size(entries,2)==2, ...
    'frameRanges must be a two-column cell array.');
for i = 1:size(entries,1)
    token = normaliseKey(entries{i,1});
    if strcmp(token, key) || strcmp(token, dataset)
        r = double(entries{i,2});
        if isvector(r), r = reshape(r,1,2); end
        assert(size(r,2)==2, 'Each frame range must be Nx2.');
        ranges = [ranges; sort(r,2)]; %#ok<AGROW>
    end
end
end

function tf = inRanges(values, ranges)
tf = false(size(values));
for i = 1:size(ranges,1)
    tf = tf | (values >= ranges(i,1) & values <= ranges(i,2));
end
end

function txt = rangesText(ranges)
if isempty(ranges), txt = "all"; return; end
parts = strings(size(ranges,1),1);
for i = 1:size(ranges,1)
    parts(i) = sprintf('[%g,%g]', ranges(i,1), ranges(i,2));
end
txt = strjoin(parts, ';');
end

function frames = absoluteFrames(R, n)
if isfield(R, 'commonFrames') && ~isempty(R.commonFrames)
    frames = nan(n,1);
    src = R.commonFrames(:);
    frames(1:min(n,numel(src))) = src(1:min(n,numel(src)));
else
    frames = (1:n).';
end
end

function animal = inferAnimal(dataset)
token = regexp(char(dataset), '^(.*?)(?:Top|Bottom|Run|Flick)', 'tokens', 'once');
if isempty(token), animal = string(dataset); else, animal = string(token{1}); end
end

function value = getScalarField(S, field, fallback)
if isfield(S, field) && isscalar(S.(field)) && isfinite(S.(field))
    value = double(S.(field));
else
    value = fallback;
end
end

function q = percentile(x, p)
x = sort(x(isfinite(x)));
if isempty(x), q = NaN; return; end
if isscalar(x), q=x; return; end
pos = 1 + (numel(x)-1) * p/100;
lo = floor(pos); hi = ceil(pos);
q = x(lo) + (pos-lo) * (x(hi)-x(lo));
end

function [group, dataset] = splitDatasetKey(key)
parts = strsplit(normaliseKey(key), '/');
assert(numel(parts)>=2, 'Dataset keys must be group/dataset (got %s).', key);
group = parts{end-1}; dataset = parts{end};
end

function list = normaliseList(list)
if ischar(list) || isstring(list), list = cellstr(string(list)); end
assert(iscell(list), 'datasetList must be a cell or string array.');
list = cellfun(@normaliseKey, list, 'UniformOutput', false);
end

function key = normaliseKey(key)
key = char(strtrim(string(key)));
key = strrep(key, '\', '/');
while contains(key, '//'), key = strrep(key, '//', '/'); end
key = strip(key, '/');
end

function T = vertcatSafe(rows)
if isempty(rows), T = table(); else, T = vertcat(rows{:}); end
end

function opts = applyDefaults(opts, defaults)
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts,names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end

function s = onoff(tf)
if tf, s='on'; else, s='off'; end
end

function s = upperFirst(s)
if ~isempty(s), s(1)=upper(s(1)); end
end
