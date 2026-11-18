function [errors, cheetahpoints] = cheetahError( ...
    params0, Rz, Rx, Ry, ...
    scaleHead, scaleTorsoLength, ...
    undistLabels, poses, intrinsicParams, LLH, ...
    frameNum, variances, params_before, camIDs, useTemporal, robustOpt, camMask)

    % robustOpt (optional): per-observation outlier handling so one inconsistent
    % camera cannot drag the pose. Defaults to plain least squares when omitted,
    % so existing call sites (e.g. the diagnostic reprojection-RMSE passes) keep
    % their original raw-residual behaviour.
    %   robustOpt.type : 'none' | 'huber' | 'tukey' | 'hard'
    %   robustOpt.c    : threshold in SIGMA units (residuals are (proj-label)/sigma)
    if nargin < 16 || isempty(robustOpt)
        robustOpt = struct('type','none','c',3);
    end
    if ~isfield(robustOpt,'type') || isempty(robustOpt.type), robustOpt.type = 'none'; end
    if ~isfield(robustOpt,'c')    || isempty(robustOpt.c),    robustOpt.c    = 3;      end

    % camMask (optional): per-frame logical over camIDs. A false entry HARD-EXCLUDES
    % that camera's observations from this frame's cost (treated exactly like an
    % invisible landmark, LLH=0), so a camera flagged as a systematic outlier by the
    % leave-one-out trust pass (cameraTrustLOO.m) cannot drag the pose. Omitted
    % or empty => all cameras used, so the raw diagnostic passes are unaffected.
    if nargin < 17 || isempty(camMask)
        camMask = true(1, numel(camIDs));
    end

    % --- Per-landmark CONSENSUS outlier rejection -----------------------------
    % robustOpt.c above is an ABSOLUTE (sigma-unit) threshold and cannot tell a
    % legitimately precise close-up camera (small residual) from a systematically
    % biased one (e.g. a mis-calibrated / sub-frame-desynced view whose error is
    % AMPLIFIED when the animal is large in it). The consensus rule judges each
    % camera RELATIVE to the other cameras viewing the SAME landmark this frame,
    % so a view that grossly disagrees with its peers is down-weighted / rejected
    % regardless of absolute pixel scale, while consistent views (close or far)
    % are kept. Smooth (redescending) so it does not destabilise lsqnonlin, and
    % only applied when >= minCams cameras see the landmark.
    %   .enable : default ON whenever a robust loss is active (type ~= 'none');
    %             the raw diagnostic passes (robustOpt absent -> 'none') keep all.
    %   .C      : reject a camera whose residual exceeds C x the per-landmark
    %             consensus scale (sigma units).
    %   .floorSigma : consensus scale floor, so well-fit frames are never trimmed.
    if ~isfield(robustOpt,'consensus') || ~isstruct(robustOpt.consensus)
        robustOpt.consensus = struct();
    end
    cc = robustOpt.consensus;
    if ~isfield(cc,'enable')     || isempty(cc.enable),     cc.enable     = ~strcmp(robustOpt.type,'none'); end
    if ~isfield(cc,'C')          || isempty(cc.C),          cc.C          = 4;  end
    if ~isfield(cc,'floorSigma') || isempty(cc.floorSigma), cc.floorSigma = 2;  end
    if ~isfield(cc,'minCams')    || isempty(cc.minCams),    cc.minCams    = 3;  end
    useConsensus = cc.enable;

    ncams = numel(camIDs);

    [~, ~, cheetahpoints] = cheetah3D( ...
        params0, Rz, Rx, Ry, ...
        scaleHead, scaleTorsoLength);

    % Order MUST match varsUsed / LLH tensor in the main script
    varnames = {'l_eye','r_eye','nose','l_ear','r_ear','spine','neck_base'};
    nparts   = numel(varnames);

    % Temporal block is 6 DOF pairs (head 6 + torso 6), independent of nparts.
    % Size the error array to accommodate whichever is larger.
    nrows = max(nparts, 6);

    error = struct;

    % ---- Pass 1: raw normalised residuals (proj-label)/sigma per (part, camera) ----
    rawErr = zeros(nparts, ncams, 2);
    gateOK = false(nparts, ncams);
    behindPenaltyPx = [1000 1000];

    for ii = 1:ncams
        camname = sprintf('cam%d', camIDs(ii));
        camOK   = camMask(ii);   % false => camera hard-excluded this frame (LOO trust)
        for i = 1:nparts
            bodypart = varnames{i};

            coord_cheetah = cheetahpoints.(bodypart)';
            pixels = world2img(coord_cheetah, invert(poses{ii}), intrinsicParams{ii});
            pix1   = undistLabels.(camname).(bodypart){frameNum};

            sigma_good = variances.(bodypart);
            sigma_mid  = sigma_good * 3;

            % Gate on the LABEL only. A labelled point whose model projection
            % leaves the image (or falls behind the camera) must still cost,
            % otherwise the fit can lower its error by pushing the model off a
            % camera's image instead of fitting it.
            if camOK && LLH(frameNum, i, ii) > 0.5 && all(isfinite(pix1)) && all(pix1 > 1)
                camPt = transformPointsForward(invert(poses{ii}), coord_cheetah);
                if camPt(3) <= 0 || ~all(isfinite(pixels))
                    pixels = pix1 + behindPenaltyPx;   % behind the camera: large constant miss
                end
                if LLH(frameNum, i, ii) < 0.85
                    rawErr(i,ii,:) = (pixels - pix1) * (1 / sigma_mid);
                else
                    rawErr(i,ii,:) = (pixels - pix1) * (1 / sigma_good);
                end
                gateOK(i,ii) = true;
            end
        end
    end

    % ---- Pass 2: per-landmark consensus weights (reject gross peer-outlier views) ----
    consW = ones(nparts, ncams);
    if useConsensus
        for i = 1:nparts
            cams = find(gateOK(i,:));
            if numel(cams) >= cc.minCams
                mags  = sqrt(rawErr(i,cams,1).^2 + rawErr(i,cams,2).^2);
                mags  = mags(:).';
                scale = max(median(mags), cc.floorSigma);
                for jj = 1:numel(cams)
                    ratio = mags(jj) / scale;
                    if ratio <= cc.C
                        consW(i,cams(jj)) = (1 - (ratio/cc.C)^2)^2;   % redescending
                    else
                        consW(i,cams(jj)) = 0;                        % reject outlier view
                    end
                end
            end
        end
    end

    % ---- Pass 3: assemble per-camera residual blocks ----
    % Per-observation absolute robust loss (huber/tukey/'none') as a backstop for
    % landmarks seen by < minCams (no consensus available), then multiply by the
    % per-landmark consensus weight.
    for ii = 1:ncams
        errVect = zeros(nrows, 2);
        for i = 1:nparts
            if gateOK(i,ii)
                err = [rawErr(i,ii,1), rawErr(i,ii,2)];
                err = robustifyResidual(err, robustOpt);
                err = err * sqrt(consW(i,ii));
                errVect(i,:) = err;
            end
        end
        error.(sprintf('camera%d', ii)) = errVect;
    end

    if useTemporal
        weightsHead  = [0.08 0.08 0.08 0.7 0.7 0.7];
        weightsTorso = [0.08 0.08 0.08 0.7 0.7 0.7];

        angle_diff = @(a, b) mod(a - b + pi, 2*pi) - pi;

        temporalHeadrpy  = angle_diff(params0(1:3), params_before(1:3));
        temporalTorsorpy = angle_diff(params0(7:9), params_before(7:9));

        temporalHeadxyz  = params0(4:6)   - params_before(4:6);
        temporalTorsoxyz = params0(10:12) - params_before(10:12);

        temporalHead  = [temporalHeadrpy  temporalHeadxyz]  ./ weightsHead;
        temporalTorso = [temporalTorsorpy temporalTorsoxyz] ./ weightsTorso;

        k = [temporalHead(:), temporalTorso(:)];  % 6x2
        error.temporal = zeros(nrows, 2);
        error.temporal(1:6, :) = k;
    else
        error.temporal = zeros(nrows, 2);
    end

    errors = zeros(nrows, 2, ncams + 1);

    for ii = 1:ncams
        errors(:,:,ii) = error.(sprintf('camera%d', ii));
    end

    errors(:,:,ncams+1) = error.temporal;
end

function err = robustifyResidual(err, ro)
    % Scale a 2-D normalised residual err = (proj-label)/sigma by sqrt(w), where
    % w is a robust weight on the residual magnitude r = |err| (in sigma units).
    % Minimising sum((sqrt(w).*err).^2) is then equivalent to minimising the
    % corresponding robust loss (IRLS embedded in lsqnonlin's Gauss-Newton step).
    if strcmp(ro.type, 'none')
        return
    end
    r = sqrt(err(1)^2 + err(2)^2);
    if r == 0
        return
    end
    switch ro.type
        case 'huber'   % gentle: quadratic within c, linear (down-weighted) beyond
            if r <= ro.c, w = 1; else, w = ro.c / r; end
        case 'tukey'   % redescending: smoothly vanishes beyond c (outlier ignored)
            if r <= ro.c, w = (1 - (r/ro.c)^2)^2; else, w = 0; end
        case 'hard'    % hard reject beyond c
            if r <= ro.c, w = 1; else, w = 0; end
        otherwise
            w = 1;
    end
    err = err * sqrt(w);
end
