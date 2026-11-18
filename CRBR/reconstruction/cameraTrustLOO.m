function [camKeptMask, camLOOResid] = cameraTrustLOO( ...
        obsCell, cameraPoses, intrinsicsForTri, posesPlot, intrinsicParams, ...
        camIDs, trustOpt)
%CAMERATRUSTLOO  Per-frame, per-camera trust via leave-one-out reprojection.
%
% Detects cameras carrying a SYSTEMATIC (coherent, proximity-amplified)
% calibration / sub-frame-sync bias, which an ordinary in-sample reprojection
% residual cannot see: when one camera is biased, the joint triangulation moves
% toward it, LOWERING its own residual and raising the peers'. The discriminating
% signal is leave-one-out (LOO): triangulate each landmark from the OTHER cameras
% and reproject into the held-out camera. A biased camera then stands out even
% when its in-sample residual looks small.
%
% INPUTS
%   obsCell          : {nframes x nBody} cell. Each entry is [] or a struct with
%                        .viewIDs (Mx1 uint/double, actual camera numbers) and
%                        .pts2D   (Mx2 undistorted pixel coords), one row per view.
%   cameraPoses      : table(ViewId, AbsolutePose) for triangulateMultiview.
%   intrinsicsForTri : intrinsics array passed to triangulateMultiview.
%   posesPlot        : {1 x ncams} camera->world rigidtform3d, in camIDs order.
%   intrinsicParams  : {1 x ncams} cameraIntrinsics, in camIDs order.
%   camIDs           : 1 x ncams actual camera numbers (column order of outputs).
%   trustOpt         : struct of options (any missing field gets the default):
%       .C_rel          (3)  reject cam k if LOOscore(k) > C_rel * median(peers)
%       .floorPx        (8)  ... AND > floorPx (so well-fit frames keep all cams)
%       .minPresentCams (4)  only run rejection when >= this many cams are present
%       .keepAtLeast    (3)  never let the kept set drop below this many cams
%       .H              (2)  hysteresis: flip kept<->rejected only after H
%                            consecutive frames (avoids frame-to-frame flicker)
%       .minCamsLOO     (3)  a landmark needs >= this many views for a LOO test
%                            (so the leave-one-out subset keeps >= 2 views)
%
% OUTPUTS
%   camKeptMask : nframes x ncams logical. false  => camera REJECTED that frame
%                 (its observations should be excluded from triangulation / fit).
%                 Cameras with no observations that frame are left true (harmless;
%                 they contribute nothing) -- inspect camLOOResid to tell them apart.
%   camLOOResid : nframes x ncams double. Per-camera robust (median) LOO pixel
%                 residual that frame; NaN where no LOO score was computable.
%
% Pure / decoupled: depends only on pointTrack, triangulateMultiview, world2img,
% invert (the same primitives multiviewRecon.m already uses).

    if nargin < 7 || isempty(trustOpt), trustOpt = struct(); end
    trustOpt = setDefault(trustOpt, 'C_rel',          3);
    trustOpt = setDefault(trustOpt, 'floorPx',        8);
    trustOpt = setDefault(trustOpt, 'minPresentCams', 4);
    trustOpt = setDefault(trustOpt, 'keepAtLeast',    3);
    trustOpt = setDefault(trustOpt, 'H',              2);
    trustOpt = setDefault(trustOpt, 'minCamsLOO',     3);

    [nframes, nBody] = size(obsCell);
    ncams = numel(camIDs);

    camKeptMask = true(nframes, ncams);
    camLOOResid = nan(nframes, ncams);

    % Hysteresis state (sticky), carried across frames.
    state   = true(1, ncams);    % true = kept, false = rejected
    badRun  = zeros(1, ncams);   % consecutive outlier frames
    goodRun = zeros(1, ncams);   % consecutive non-outlier frames

    for f = 1:nframes
        % ---- collect per-camera LOO residuals across this frame's landmarks ----
        looByCam = cell(1, ncams);
        present  = false(1, ncams);

        for b = 1:nBody
            o = obsCell{f,b};
            if isempty(o) || ~isstruct(o), continue; end
            vID = double(o.viewIDs(:));
            p2D = o.pts2D;
            nv  = numel(vID);

            % every viewing camera counts as "present" this frame
            for r = 1:nv
                ci = find(camIDs == vID(r), 1);
                if ~isempty(ci), present(ci) = true; end
            end

            if nv < trustOpt.minCamsLOO, continue; end

            % leave-one-out: triangulate from the OTHER views, reproject into r
            for r = 1:nv
                keepR        = true(nv,1);
                keepR(r)     = false;
                othID        = vID(keepR);
                othP         = p2D(keepR,:);
                if numel(othID) < 2, continue; end

                track = pointTrack(uint32(othID(:).'), othP);
                try
                    [X, ~, validIdx] = triangulateMultiview(track, cameraPoses, intrinsicsForTri);
                catch
                    continue
                end
                if isempty(X) || ~all(isfinite(X(1,:))) || isempty(validIdx) || ~validIdx(1)
                    continue
                end

                ci = find(camIDs == vID(r), 1);
                if isempty(ci), continue; end
                uv_proj = world2img(X(1,:), invert(posesPlot{ci}), intrinsicParams{ci});
                resid   = hypot(uv_proj(1) - p2D(r,1), uv_proj(2) - p2D(r,2));
                looByCam{ci}(end+1) = resid; %#ok<AGROW>
            end
        end

        % ---- robust per-camera LOO score (median over its landmarks) ----
        looScore = nan(1, ncams);
        for ci = 1:ncams
            if ~isempty(looByCam{ci})
                looScore(ci) = median(looByCam{ci});
            end
        end
        camLOOResid(f,:) = looScore;

        % ---- consensus outlier flags (relative threshold + floor + guards) ----
        isOut = consensusOutliers(looScore, present, trustOpt);

        % ---- symmetric hysteresis: require H consecutive frames to flip ----
        for ci = 1:ncams
            if ~present(ci)
                badRun(ci) = 0; goodRun(ci) = 0; state(ci) = true;
            elseif isOut(ci)
                badRun(ci) = badRun(ci) + 1; goodRun(ci) = 0;
                if badRun(ci) >= trustOpt.H, state(ci) = false; end
            else
                goodRun(ci) = goodRun(ci) + 1; badRun(ci) = 0;
                if goodRun(ci) >= trustOpt.H, state(ci) = true; end
            end
        end

        % ---- safety cap: keep >= keepAtLeast present cams (drop worst first) ----
        rej    = present & ~state;
        maxRej = max(0, sum(present) - trustOpt.keepAtLeast);
        if sum(rej) > maxRej
            cand = find(rej);
            [~, ord] = sort(looScore(cand), 'descend');   % worst LOO first
            forcedKeep = cand(ord(maxRej+1:end));
            state(forcedKeep) = true;                      % re-admit the least-bad
        end

        keep = true(1, ncams);
        keep(present & ~state) = false;
        camKeptMask(f,:) = keep;
    end
end

% ------------------------------------------------------------------------------
function isOut = consensusOutliers(looScore, present, t)
    ncams = numel(looScore);
    isOut = false(1, ncams);

    if sum(present) < t.minPresentCams, return; end

    scored = find(present & isfinite(looScore));
    if numel(scored) < 3, return; end   % need enough peers for a consensus

    ratio = zeros(1, ncams);
    for jj = scored
        others = looScore(setdiff(scored, jj));
        med    = median(others);
        thr    = max(t.C_rel * med, t.floorPx);
        if looScore(jj) > thr
            isOut(jj) = true;
            ratio(jj) = looScore(jj) / max(med, eps);
        end
    end

    % never flag more than (present - keepAtLeast); keep the worst offenders
    maxRej = max(0, sum(present) - t.keepAtLeast);
    if sum(isOut) > maxRej
        cand = find(isOut);
        [~, ord] = sort(ratio(cand), 'descend');
        isOut = false(1, ncams);
        isOut(cand(ord(1:maxRej))) = true;
    end
end

% ------------------------------------------------------------------------------
function s = setDefault(s, field, val)
    if ~isfield(s, field) || isempty(s.(field)), s.(field) = val; end
end
