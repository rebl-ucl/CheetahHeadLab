function [lead,rho,point]=spatialForwardExtension(origin,direction,endpoint,u,endS,currentS)
% Intersect two forward rays: head direction and the target endpoint extension.
% No finite extrapolation cap. Parallel or backward intersections stay missing.
lead=NaN;rho=NaN;point=[NaN NaN NaN];
if any(~isfinite([origin direction endpoint u endS currentS])) || norm(direction)<1e-12,return;end
d=direction/norm(direction);q=endpoint(1:2)-origin;
den=d(1)*u(2)-d(2)*u(1);
if abs(den)<=1e-10,return;end
r=(q(1)*u(2)-q(2)*u(1))/den;
a=(q(1)*d(2)-q(2)*d(1))/den;
if r<0 || a<0,return;end
rho=r;point=[endpoint(1:2)+a*u endpoint(3)];lead=endS+a-currentS;
end
