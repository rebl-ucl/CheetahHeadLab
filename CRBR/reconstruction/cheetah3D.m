function [R1_2, R1_3, cheetahpoints] = cheetah3D( ...
    params0, Rz, Rx, Ry, ...
    scaleHead, scaleTorsoLength)
% Full Head Trajectory Estimation variant: shoulders removed from kinematic
% model; left/right ears added as head landmarks for reconstruction.
% Torso width is not a free scale parameter (fixed at template proportions).
% Returns R1_2 (head rotation, local->world), R1_3 (torso rotation, local->world),
% and cheetahpoints (the world-frame landmark positions).

    params0 = num2cell(params0);
    [theta_x, theta_y, theta_z, x, y, z, ...
     theta_x_2, theta_y_2, theta_z_2, x2, y2, z2] = deal(params0{:});

    cheetahpoints = struct;

    %% ---------------- Rotations ----------------
    R1_2 = Rz(theta_z)   * Rx(theta_x)   * Ry(theta_y);     % head
    R1_3 = Rz(theta_z_2) * Rx(theta_x_2) * Ry(theta_y_2);   % torso

    %% ---------------- Origins ----------------
    phead  = [x;  y;  z];
    ptorso = [x2; y2; z2];

    %% ---------------- Head geometry ----------------
    % +Y = left, -Y = right
    pleye = phead + R1_2 * [0;  scaleHead*0.03; 0];
    preye = phead + R1_2 * [0; -scaleHead*0.03; 0];

    % Ears: backward (-X) and above (+Z) the centre of eyes
    plear = phead + R1_2 * [scaleHead*-0.0725;  scaleHead* 0.0608; scaleHead* 0.0441];
    prear = phead + R1_2 * [scaleHead*-0.0725; scaleHead*-0.0608; scaleHead* 0.0441];

    pnose = phead + R1_2 * [scaleHead*0.055; 0; -scaleHead*0.055];

    %% ---------------- Torso geometry (no shoulders) ----------------
    pneck_base = ptorso + R1_3 * [scaleTorsoLength*0.04;  0; scaleTorsoLength*0.10];
    pspine_mid = ptorso + R1_3 * [-scaleTorsoLength*0.33; 0; scaleTorsoLength*0.10];

    %% ---------------- Store ----------------
    cheetahpoints.l_eye     = pleye;
    cheetahpoints.r_eye     = preye;
    cheetahpoints.nose      = pnose;
    cheetahpoints.l_ear     = plear;
    cheetahpoints.r_ear     = prear;
    cheetahpoints.spine     = pspine_mid;
    cheetahpoints.neck_base = pneck_base;
end
