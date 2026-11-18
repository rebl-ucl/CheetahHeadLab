function [y, nSpike] = mcFilterAngleSeries(x, mask, K, nSigma, sgOrder, sgWin, isUnsigned)
%MCFILTERANGLESERIES Hampel despike + Savitzky-Golay smoothing for MC angles.
%   Signed component errors are re-wrapped to +/-pi. Unsigned total angles
%   are clamped to [0,pi] after smoothing so filter overshoot cannot make them
%   negative. NaN/excluded frames remain untouched.

if nargin < 7
    isUnsigned = false;
end

y = x;
nSpike = 0;
idx = find(mask(:));
if numel(idx) < 3
    return
end
seg = x(idx);

if numel(seg) >= 2*K + 1
    [seg, outliers] = hampel(seg, K, nSigma);
    nSpike = sum(outliers);
end

w = sgWin;
if mod(w, 2) == 0
    w = w - 1;
end
if numel(seg) >= w && w > sgOrder
    seg = sgolayfilt(seg, sgOrder, w);
end

if isUnsigned
    y(idx) = min(max(seg, 0), pi);
else
    y(idx) = wrapToPi(seg);
end
end
