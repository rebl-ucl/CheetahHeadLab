function estimateScales(parent_folder, opts)
% estimateScales  CRBR head and torso model scaling from triangulation.
%
%   estimateScales(parent_folder, opts)
%
% Uses the multiview triangulation (triangulatedPoints*_multiview.mat) to
% estimate the head and torso-length scales, saves them to
% <parent_folder>/cheetah_scales.mat, and (optionally) draws the template
% model (scale 1) against the scaled model. Validation is pose-free: it
% compares median measured triangulated distances against template and scaled
% distances. Does not use camera IDs.
%
% opts fields (all optional):
%   triSource("raw")   which triangulation source to read ("raw"/"med"/"lpf")
%   showFigures(false) draw the head/torso scaling comparison figure
%
% This is consistent with the current cheetah3D geometry:
%   Head:  l_eye=[0;+0.03;0]  r_eye=[0;-0.03;0]  nose=[0.055;0;-0.055]
%   Torso: neck_base=[0.04;0;0.10]  spine=[-0.33;0;0.10]
%
% Called with no arguments it runs a default dataset (standalone debugging).

if nargin < 1 || isempty(parent_folder)
    parent_folder = fullfile(pwd, 'data', 'pitch_tracking', 'phantomTopRun13_20170902');
end
if nargin < 2 || isempty(opts)
    opts = struct();
end

triSource   = getOpt(opts, 'triSource', "raw");   % "raw", "med", "lpf"
showFigures = getOpt(opts, 'showFigures', false);

%% -------------------- Rotation helpers --------------------
Rx = @(theta_x) [1 0 0;
                 0 cos(theta_x) sin(theta_x);
                 0 -sin(theta_x) cos(theta_x)];

Ry = @(theta_y) [cos(theta_y) 0 sin(theta_y);
                 0 1 0;
                 -sin(theta_y) 0 cos(theta_y)];

Rz = @(theta_z) [cos(theta_z) sin(theta_z) 0;
                 -sin(theta_z) cos(theta_z) 0;
                 0 0 1];

%% -------------------- Load triangulated points --------------------
minViewsForScale = 2;   % minimum cameras that must agree on a point for it to count

switch lower(triSource)
    case "raw"
        S = load(fullfile(parent_folder, "triangulatedPointsRaw_multiview.mat"), ...
            "pointsRaw", "numViewsUsed", "varsUsed");
        pointsStruct = S.pointsRaw;
    case "med"
        S = load(fullfile(parent_folder, "triangulatedPointsMed_multiview.mat"), ...
            "pointsMed", "numViewsUsed", "varsUsed");
        pointsStruct = S.pointsMed;
    case "lpf"
        S = load(fullfile(parent_folder, "triangulatedPointsLPF_multiview.mat"), ...
            "pointsLPF", "numViewsUsed", "varsUsed");
        pointsStruct = S.pointsLPF;
    otherwise
        error('triSource must be "raw", "med", or "lpf"');
end

numViewsUsed = S.numViewsUsed;   % [nframes x nBody]
varsInFile   = S.varsUsed;

nframes = size(pointsStruct.nose, 1);

% Helper: per-frame camera count for a named body part
getViews = @(bp) numViewsUsed(:, strcmp(varsInFile, bp));

%% -------------------- Compute triangulated distances across frames --------------------
leye  = pointsStruct.l_eye;
reye  = pointsStruct.r_eye;
nose  = pointsStruct.nose;
neck  = pointsStruct.neck_base;
spine = pointsStruct.spine;

coe = (leye + reye) / 2;

d_eye       = vecnorm(leye - reye, 2, 2);     % eye width
d_nose      = vecnorm(nose - coe, 2, 2);      % coe -> nose
d_neckspine = vecnorm(neck - spine, 2, 2);    % neck -> spine

% Each scale uses only the body parts it actually needs,
% gated by a minimum camera count to reject poorly-triangulated frames
valid_head        = isfinite(d_eye) & isfinite(d_nose) & ...
                    getViews('l_eye') >= minViewsForScale & ...
                    getViews('r_eye') >= minViewsForScale & ...
                    getViews('nose')  >= minViewsForScale;
valid_torsoLength = isfinite(d_neckspine) & ...
                    getViews('neck_base') >= minViewsForScale & ...
                    getViews('spine')     >= minViewsForScale;
assert(any(valid_head),        "No valid frames for head scale (try lowering minViewsForScale).");
assert(any(valid_torsoLength), "No valid frames for torso length (try lowering minViewsForScale).");

fprintf('Frames passing validity gate: head=%d  torsoLen=%d  (of %d total)\n', ...    
    sum(valid_head), sum(valid_torsoLength), nframes);    

% robust medians
d_eye_med       = median(d_eye(valid_head),               'omitnan');
d_nose_med      = median(d_nose(valid_head),              'omitnan');
d_neckspine_med = median(d_neckspine(valid_torsoLength),  'omitnan');

%% -------------------- Template distances from unscaled model --------------------
% These MUST match cheetah3D(..., scaleHead=1, scaleTorsoLength=1)

% ----- Head template -----
template_leye = [0;  0.03; 0];
template_reye = [0; -0.03; 0];
template_coe  = [0;  0;    0];
template_nose = [0.055; 0; -0.055];

template_eye       = norm(template_leye - template_reye);
template_nose_dist = norm(template_nose - template_coe);

% ----- Torso template -----
template_neck  = [ 0.04;  0.00; 0.10];
template_spine = [-0.33;  0.00; 0.10];

template_neckspine = norm(template_neck - template_spine);

% Diagnostic: distribution of measured distances vs template values
fprintf('\nd_eye      distribution: min=%.3f  p25=%.3f  median=%.3f  p75=%.3f  max=%.3f  (template=%.3f)\n', ...    
    min(d_eye(valid_head)),         prctile(d_eye(valid_head), 25), ...  
    median(d_eye(valid_head)),      prctile(d_eye(valid_head), 75), ...  
    max(d_eye(valid_head)),         template_eye);    
fprintf('d_neckspine distribution: min=%.3f  p25=%.3f  median=%.3f  p75=%.3f  max=%.3f  (template=%.3f)\n\n', ...    
    min(d_neckspine(valid_torsoLength)),    prctile(d_neckspine(valid_torsoLength), 25), ...  
    median(d_neckspine(valid_torsoLength)), prctile(d_neckspine(valid_torsoLength), 75), ...  
    max(d_neckspine(valid_torsoLength)),    template_neckspine);    

%% -------------------- Estimate scales --------------------
s_eye       = d_eye_med       / template_eye;
s_nose      = d_nose_med      / template_nose_dist;
s_neckspine = d_neckspine_med / template_neckspine;

% Head scale: always the nose-based estimate. Eye width is the smallest-
% baseline landmark pair (template 0.06 m), so it is hypersensitive to label
% noise, weak camera coverage or an L/R swap, and across the datasets it runs
% about 20% above s_nose (the labelled eye spacing is wider than the
% template). Averaging it in only when the two agree within a threshold made
% the head scale jump by ~10% depending on which side of the threshold a run
% fell, so s_eye is now a diagnostic only.
headScaleDivergence = abs(s_eye - s_nose) / mean([s_eye, s_nose]);
scaleHead     = s_nose;
headScaleNote = sprintf('  (= s_nose; s_eye differs by %.0f%%, diagnostic only)', ...
    headScaleDivergence * 100);
scaleTorsoLength = s_neckspine;

fprintf('\nEstimated scales from triangulation:\n');
fprintf('s_eye            = %.4f\n', s_eye);
fprintf('s_nose           = %.4f\n', s_nose);
fprintf('s_neckspine      = %.4f\n', s_neckspine);
fprintf('scaleHead        = %.4f%s\n', scaleHead, headScaleNote);
fprintf('scaleTorsoLength = %.4f\n\n', scaleTorsoLength);

%% -------------------- Save results --------------------
save(fullfile(parent_folder, "cheetah_scales.mat"), ...
    "scaleHead", "scaleTorsoLength");

%% -------------------- Build nominal-pose models (no recon pose needed) --------------------
% Visualise the SCALING effect on geometry, not the world pose, so use a
% default nominal pose: identity rotations, origin at zero. Head and torso
% are drawn in their OWN local frames (two subplots), so no reconstructed pose
% is required to place them.
nominalParams = zeros(1, 12);   % [head rot(3) xyz(3), torso rot(3) xyz(3)]

[~, ~, pts_before] = cheetah3D(nominalParams, Rz, Rx, Ry, 1, 1);
[~, ~, pts_after]  = cheetah3D(nominalParams, Rz, Rx, Ry, ...
                               scaleHead, scaleTorsoLength);

headParts  = {'l_eye','r_eye','nose','l_ear','r_ear'};
torsoParts = {'neck_base','spine'};

%% -------------------- Plot head and torso in separate local frames --------------------
if showFigures
figure('Color','w');
set(gcf, 'Position', [100 100 1400 600]);

% ----- Head (head-local frame) -----
subplot(1,2,1); hold on; grid on; axis equal; view(3);
plotSkeletonHead3D(pts_before, [1 0 0], 2.5, '--');
plotPointSet3D(pts_before, headParts, [1 0 0], 50, false);
plotSkeletonHead3D(pts_after,  [0 0.6 0], 2.5, '-');
plotPointSet3D(pts_after,  headParts, [0 0.6 0], 50, false);
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Head (local frame, nominal pose)\nscaleHead = %.3f', scaleHead));
legend({'Template (scale 1)', 'Scaled'}, 'Location', 'best');

% ----- Torso (torso-local frame) -----
subplot(1,2,2); hold on; grid on; axis equal; view(3);
plotSkeletonTorso3D(pts_before, [1 0 0], 2.5, '--');
plotPointSet3D(pts_before, torsoParts, [1 0 0], 50, false);
plotSkeletonTorso3D(pts_after,  [0 0.6 0], 2.5, '-');
plotPointSet3D(pts_after,  torsoParts, [0 0.6 0], 50, false);
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Torso (local frame, nominal pose)\nscaleTorsoLength = %.3f', scaleTorsoLength));
legend({'Template (scale 1)', 'Scaled'}, 'Location', 'best');
end

%% -------------------- Pose-free validation of scaling --------------------
% Compare the median measured triangulated distances against the template
% (scale 1) and the scaled-model distances. No pose required: scaling acts
% on scalar distances, so a good scale brings "scaled" close to "measured".
fprintf('Scaling check (median measured vs template vs scaled model):\n');
fprintf('%-12s  %9s  %9s  %9s\n', 'distance', 'measured', 'template', 'scaled');
fprintf('%-12s  %9.4f  %9.4f  %9.4f\n', 'eye width',  d_eye_med,       template_eye,       template_eye*scaleHead);
fprintf('%-12s  %9.4f  %9.4f  %9.4f\n', 'nose',       d_nose_med,      template_nose_dist, template_nose_dist*scaleHead);
fprintf('%-12s  %9.4f  %9.4f  %9.4f\n', 'neck-spine', d_neckspine_med, template_neckspine, template_neckspine*scaleTorsoLength);

end

%% -------------------- Local functions --------------------
function v = getOpt(opts, name, default)
%GETOPT  Return opts.(name) if the field exists, else default.
    if isstruct(opts) && isfield(opts, name)
        v = opts.(name);
    else
        v = default;
    end
end

function plotPointSet3D(ptsStruct, varsUsed, colorVal, markerSize, filledFlag)
    first = true;
    for i = 1:numel(varsUsed)
        nm = varsUsed{i};
        if isfield(ptsStruct, nm)
            p = ptsStruct.(nm)(:).';
            if all(isfinite(p))
                if filledFlag
                    if first
                        scatter3(p(1), p(2), p(3), markerSize, colorVal, ...
                            'filled', 'DisplayName', 'Triangulated points');
                        first = false;
                    else
                        scatter3(p(1), p(2), p(3), markerSize, colorVal, ...
                            'filled', 'HandleVisibility','off');
                    end
                else
                    scatter3(p(1), p(2), p(3), markerSize, colorVal, ...
                        'LineWidth', 1.5, 'HandleVisibility','off');
                end
                text(p(1), p(2), p(3), ['  ' nm], 'Color', colorVal, 'FontSize', 10);
            end
        end
    end
end

function plotSkeletonHead3D(ptsStruct, colorVal, lw, lineStyle)
    drawSeg(ptsStruct, 'l_eye', 'r_eye', colorVal, lw, lineStyle);
    drawMidSeg(ptsStruct, 'l_eye', 'r_eye', 'nose', colorVal, lw, lineStyle);
    drawSeg(ptsStruct, 'l_ear', 'r_ear', colorVal, lw, lineStyle);
    drawSeg(ptsStruct, 'l_ear', 'l_eye', colorVal, lw, lineStyle);
    drawSeg(ptsStruct, 'r_ear', 'r_eye', colorVal, lw, lineStyle);
    addSkeletonLegendEntry(colorVal, lw, lineStyle);
end

function plotSkeletonTorso3D(ptsStruct, colorVal, lw, lineStyle)
    drawSeg(ptsStruct, 'neck_base', 'spine', colorVal, lw, lineStyle);
    addSkeletonLegendEntry(colorVal, lw, lineStyle);
end

function addSkeletonLegendEntry(colorVal, lw, lineStyle)
    if isequal(colorVal, [1 0 0])
        plot3(nan, nan, nan, 'Color', colorVal, 'LineStyle', lineStyle, ...
              'LineWidth', lw, 'DisplayName', 'Template (scale 1)');
    elseif isequal(colorVal, [0 0.6 0])
        plot3(nan, nan, nan, 'Color', colorVal, 'LineStyle', lineStyle, ...
              'LineWidth', lw, 'DisplayName', 'Scaled');
    end
end

function drawSeg(ptsStruct, a, b, colorVal, lw, lineStyle)
    if isfield(ptsStruct, a) && isfield(ptsStruct, b)
        pa = ptsStruct.(a)(:).';
        pb = ptsStruct.(b)(:).';
        if all(isfinite(pa)) && all(isfinite(pb))
            plot3([pa(1) pb(1)], [pa(2) pb(2)], [pa(3) pb(3)], ...
                'Color', colorVal, 'LineWidth', lw, ...
                'LineStyle', lineStyle, 'HandleVisibility','off');
        end
    end
end

function drawMidSeg(ptsStruct, a, b, c, colorVal, lw, lineStyle)
    if isfield(ptsStruct, a) && isfield(ptsStruct, b) && isfield(ptsStruct, c)
        pa = ptsStruct.(a)(:).';
        pb = ptsStruct.(b)(:).';
        pc = ptsStruct.(c)(:).';
        if all(isfinite(pa)) && all(isfinite(pb)) && all(isfinite(pc))
            pm = (pa + pb)/2;
            plot3([pm(1) pc(1)], [pm(2) pc(2)], [pm(3) pc(3)], ...
                'Color', colorVal, 'LineWidth', lw, ...
                'LineStyle', lineStyle, 'HandleVisibility','off');
        end
    end
end
