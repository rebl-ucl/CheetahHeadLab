function lurePos = mcLureTriangulate(lab2D, vis, camPosesTbl, intrinArr, camIDs)
% MCLURETRIANGULATE  Re-triangulate the lure 3D position from (possibly
% perturbed) raw 2D labels for ONE frame, mirroring the production lure
% pipeline (lureReconTriangulation.m:235-265).
%
%   lurePos = mcLureTriangulate(lab2D, vis, camPosesTbl, intrinArr, camIDs)
%
% Inputs:
%   lab2D       [ncams x 2] RAW (distorted) lure pixels per camera (NaN if absent)
%   vis         [ncams x 1] logical visibility per camera
%   camPosesTbl table(ViewId uint32, AbsolutePose rigidtform3d), ALL cameras
%   intrinArr   [1 x ncams] cameraIntrinsics, aligned with camPosesTbl rows/camIDs
%   camIDs      [1 x ncams] physical camera IDs (same order as intrinArr)
%
% Output:
%   lurePos     1x3 triangulated world point, or [NaN NaN NaN] if <2 visible
%               views or triangulation fails/returns invalid.
%
% Each visible camera's raw label is undistorted with its own intrinsics
% (undistortPoints), exactly as production, then triangulateMultiview matches
% the pointTrack's ViewIds against camPosesTbl.ViewId.

    ncams = numel(camIDs);
    viewIDs = zeros(0, 1, 'uint32');
    pts2D   = zeros(0, 2);

    for ii = 1:ncams
        if ~vis(ii), continue; end
        uv = lab2D(ii, :);
        if ~all(isfinite(uv)), continue; end
        uvUD = undistortPoints(uv, intrinArr(ii));   % raw -> undistorted (production)
        viewIDs(end+1, 1) = uint32(camIDs(ii)); %#ok<AGROW>
        pts2D(end+1, :)   = uvUD;                %#ok<AGROW>
    end

    if numel(viewIDs) < 2
        lurePos = [NaN NaN NaN];
        return
    end

    track = pointTrack(viewIDs.', pts2D);
    try
        [X, ~, validIdx] = triangulateMultiview(track, camPosesTbl, intrinArr);
    catch
        lurePos = [NaN NaN NaN];
        return
    end

    if isempty(X) || ~all(isfinite(X(1, :))) || isempty(validIdx) || ~validIdx(1)
        lurePos = [NaN NaN NaN];
        return
    end
    lurePos = X(1, :);
end
