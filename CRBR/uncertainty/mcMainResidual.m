function r = mcMainResidual(params, Rz, Rx, Ry, off, ...
        labelsXY, llhMat, sigmas, invPoses, intrin, ncams, nparts, ...
        partNames, robustOpt, camMask)
% MCMAINRESIDUAL  Single-frame reprojection residual for the main CRBR analysis,
% with the model geometry given by explicit local offsets `off`.
%
% Reproduces CRBR/reconstruction/cheetahError EXACTLY (the production per-frame estimator):
% three passes -- (1) sigma-normalised reprojection residual to UNDISTORTED
% pixels with the label-only LLH gate (behind-camera points cost a 1000 px miss)
% + camMask hard-exclusion, (2) per-landmark consensus
% outlier weights, (3) assembly with per-observation robust loss x consensus
% weight. No temporal term (per-frame method). Geometry comes from
% mcGeomFromOffsets so the kinematic lengths can be perturbed per MC sample.
%
% camMask (1 x ncams logical): a false entry HARD-EXCLUDES that camera this
% frame (leave-one-out trust pass), exactly like an invisible landmark (LLH=0).
% robustOpt.consensus (struct .enable/.C/.floorSigma/.minCams): per-landmark
% consensus rejection, defaults matching cheetahError (ON whenever robust type
% ~= 'none').

    if nargin < 15 || isempty(camMask)
        camMask = true(1, ncams);
    end

    % Consensus config (mirror cheetahError.m defaults)
    if ~isfield(robustOpt, 'consensus') || ~isstruct(robustOpt.consensus)
        robustOpt.consensus = struct();
    end
    cc = robustOpt.consensus;
    if ~isfield(cc,'enable')     || isempty(cc.enable),     cc.enable     = ~strcmp(robustOpt.type,'none'); end
    if ~isfield(cc,'C')          || isempty(cc.C),          cc.C          = 4;  end
    if ~isfield(cc,'floorSigma') || isempty(cc.floorSigma), cc.floorSigma = 2;  end
    if ~isfield(cc,'minCams')    || isempty(cc.minCams),    cc.minCams    = 3;  end
    useConsensus = cc.enable;

    cp = mcGeomFromOffsets(params, Rz, Rx, Ry, off);

    % ---- Pass 1: raw normalised residuals (proj-label)/sigma per (part, camera) ----
    % All gated landmarks of a camera are projected in ONE world2img call: the
    % per-call overhead, not the arithmetic, dominates the cost. Residuals are
    % bit-identical to projecting the landmarks one at a time.
    rawErr = zeros(nparts, ncams, 2);
    gateOK = false(nparts, ncams);
    behindPenaltyPx = [1000 1000];
    Pall = zeros(nparts, 3);
    for p = 1:nparts
        Pall(p,:) = cp.(partNames{p}).';
    end
    for ii = 1:ncams
        if ~camMask(ii), continue; end   % camera hard-excluded this frame (LOO trust)
        labs = reshape(labelsXY(:, ii, :), nparts, 2);
        L = llhMat(:, ii);
        % Gate on the LABEL only (as cheetahError). A labelled point whose
        % model projection leaves the image (or falls behind the camera)
        % must still cost, otherwise the fit can lower its error by pushing
        % the model off a camera's image instead of fitting it.
        g = find(L > 0.5 & all(isfinite(labs), 2) & all(labs > 1, 2));
        if isempty(g), continue; end
        uv = world2img(Pall(g,:), invPoses{ii}, intrin{ii});
        camPt = transformPointsForward(invPoses{ii}, Pall(g,:));
        for j = 1:numel(g)
            p = g(j);
            lab = labs(p,:);
            u = uv(j,:);
            if camPt(j,3) <= 0 || ~all(isfinite(u))
                u = lab + behindPenaltyPx;   % behind the camera: large constant miss
            end
            if L(p) < 0.85
                s = sigmas(p) * 3;
            else
                s = sigmas(p);
            end
            rawErr(p,ii,:) = (u - lab) / s;
            gateOK(p,ii)   = true;
        end
    end

    % ---- Pass 2: per-landmark consensus weights (reject gross peer-outlier views) ----
    consW = ones(nparts, ncams);
    if useConsensus
        for p = 1:nparts
            cams = find(gateOK(p,:));
            if numel(cams) >= cc.minCams
                mags  = sqrt(rawErr(p,cams,1).^2 + rawErr(p,cams,2).^2);
                mags  = mags(:).';
                scale = max(median(mags), cc.floorSigma);
                for jj = 1:numel(cams)
                    ratio = mags(jj) / scale;
                    if ratio <= cc.C
                        consW(p,cams(jj)) = (1 - (ratio/cc.C)^2)^2;   % redescending
                    else
                        consW(p,cams(jj)) = 0;                        % reject outlier view
                    end
                end
            end
        end
    end

    % ---- Pass 3: assemble stacked residual (robust loss x consensus weight) ----
    r = zeros(2 * ncams * nparts, 1);
    n = 0;
    for ii = 1:ncams
        for p = 1:nparts
            n = n + 1;
            d = [0 0];
            if gateOK(p,ii)
                d = [rawErr(p,ii,1), rawErr(p,ii,2)];
                d = robustify(d, robustOpt);
                d = d * sqrt(consW(p,ii));
            end
            r(2*n-1 : 2*n) = d;
        end
    end
end

function err = robustify(err, ro)
    if strcmp(ro.type, 'none')
        return
    end
    rmag = sqrt(err(1)^2 + err(2)^2);
    if rmag == 0
        return
    end
    switch ro.type
        case 'huber'
            if rmag <= ro.c, w = 1; else, w = ro.c / rmag; end
        case 'tukey'
            if rmag <= ro.c, w = (1 - (rmag/ro.c)^2)^2; else, w = 0; end
        case 'hard'
            if rmag <= ro.c, w = 1; else, w = 0; end
        otherwise
            w = 1;
    end
    err = err * sqrt(w);
end
