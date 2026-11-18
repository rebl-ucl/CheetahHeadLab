function ang = mcMainAngles(params, lurePos, Rz, Rx, Ry)
% MCMAINANGLES  Head-stabilisation angles for the main analysis.
%
%   ang = mcMainAngles(params, lurePos, Rz, Rx, Ry)
%
% Returns (radians):
%   [alpha beta gamma delta  headLureAngle bodyLureAngle]
%
%   gamma/delta = head/body YAW error  (azimuth gaze - azimuth to-lure)
%   alpha/beta  = head/body PITCH error (elevation gaze - elevation to-lure)
%                 -- the four target-relative orientation errors.
%   headLureAngle = the TOTAL angle between the reconstructed head-forward unit
%                 vector and the unit vector from the head (coe) to the lure.
%                 This is the single "off-target" angle (combines yaw+pitch);
%                 a sanity check that should track sqrt(gamma^2+alpha^2)-ish.
%   bodyLureAngle = same for the torso.
%
% coe = params(4:6), torso_origin = params(10:12) (verified vs cheetah3D_HS_LA).
% Non-finite lure -> all NaN.

    R1_2 = Rz(params(3)) * Rx(params(1)) * Ry(params(2));
    R1_3 = Rz(params(9)) * Rx(params(7)) * Ry(params(8));

    hg = R1_2 * [1; 0; 0];   % head-forward unit vector (R orthonormal)
    tg = R1_3 * [1; 0; 0];   % torso-forward unit vector

    coe = params(4:6).';
    to  = params(10:12).';
    L   = lurePos(:);

    if all(isfinite(L))
        hl = L - coe;            % head -> lure
        tl = L - to;             % torso -> lure
        [aHG, eHG] = cart2sph(hg(1), hg(2), hg(3));
        [aTG, eTG] = cart2sph(tg(1), tg(2), tg(3));
        [aHL, eHL] = cart2sph(hl(1), hl(2), hl(3));
        [aTL, eTL] = cart2sph(tl(1), tl(2), tl(3));

        gamma = wrapToPi(aHG - aHL);
        delta = wrapToPi(aTG - aTL);
        alpha = wrapToPi(eHG - eHL);
        beta  = wrapToPi(eTG - eTL);

        uh = hl / norm(hl);
        ut = tl / norm(tl);
        headLureAngle = atan2(norm(cross(hg, uh)), dot(hg, uh));
        bodyLureAngle = atan2(norm(cross(tg, ut)), dot(tg, ut));
    else
        gamma = NaN; delta = NaN; alpha = NaN; beta = NaN;
        headLureAngle = NaN; bodyLureAngle = NaN;
    end

    ang = [alpha beta gamma delta headLureAngle bodyLureAngle];
end
