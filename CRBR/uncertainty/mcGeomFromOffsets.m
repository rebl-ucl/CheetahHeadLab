function [cp, R1_2, R1_3] = mcGeomFromOffsets(params, Rz, Rx, Ry, off)
% MCGEOMFROMOFFSETS  CRBR landmark world positions from explicit local offsets.
%
%   [cp, R1_2, R1_3] = mcGeomFromOffsets(params, Rz, Rx, Ry, off)
%
% Clean geometry for the main analysis: every head/torso landmark is
%   head:   p = phead  + R1_2 * off.<name>
%   torso:  p = ptorso + R1_3 * off.<name>
% where off.<name> are the local offset vectors IN METRES (already at the
% reconstruction size — no separate scale factor). This is the same geometry
% as CRBR/reconstruction/cheetah3D, just parameterised directly by the lengths
% so the Monte-Carlo can perturb each length independently.
%
% off fields (each a 3-vector): l_eye r_eye nose l_ear r_ear neck_base spine
% params is the 12-DOF pose [thx thy thz x y z thx2 thy2 thz2 x2 y2 z2].

    p = num2cell(params);
    [theta_x, theta_y, theta_z, x, y, z, ...
     theta_x_2, theta_y_2, theta_z_2, x2, y2, z2] = deal(p{:});

    R1_2 = Rz(theta_z)   * Rx(theta_x)   * Ry(theta_y);
    R1_3 = Rz(theta_z_2) * Rx(theta_x_2) * Ry(theta_y_2);

    phead  = [x;  y;  z];
    ptorso = [x2; y2; z2];

    cp = struct();
    cp.l_eye     = phead  + R1_2 * off.l_eye(:);
    cp.r_eye     = phead  + R1_2 * off.r_eye(:);
    cp.nose      = phead  + R1_2 * off.nose(:);
    cp.l_ear     = phead  + R1_2 * off.l_ear(:);
    cp.r_ear     = phead  + R1_2 * off.r_ear(:);
    cp.neck_base = ptorso + R1_3 * off.neck_base(:);
    cp.spine     = ptorso + R1_3 * off.spine(:);
end
