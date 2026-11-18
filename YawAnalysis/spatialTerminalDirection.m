function u=spatialTerminalDirection(frame,xyz,fs,windowSec)
% Horizontal continuation direction from a trailing least-squares line fit.
% Anchor the continuation at the final observed point, at constant height.
u=[NaN NaN];
valid=isfinite(frame(:))&all(isfinite(xyz(:,1:2)),2);
frame=double(frame(valid));xy=double(xyz(valid,1:2));
if numel(frame)<3,return;end
keep=frame>=frame(end)-windowSec*fs;
if nnz(keep)<3,keep(end-2:end)=true;end
t=(frame(keep)-frame(end))/fs;
B=[ones(numel(t),1),t]\xy(keep,:);v=B(2,:);
if norm(v)>1e-9,u=v/norm(v);end
end
