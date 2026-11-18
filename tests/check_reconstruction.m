function check_reconstruction(workFolder)
%CHECK_RECONSTRUCTION Check complete coverage and agreement with saved poses.
% check_reconstruction('outputs/julesBottomFlick2_20171209')
% Thresholds are regression tolerances, not estimates of biological uncertainty.
S=jsondecode(fileread(fullfile(workFolder,'reproduction_summary.json')));
assert(S.identical_frame_coverage,'Frame coverage differs from the production reference.');
assert(S.max_head_axis_difference_deg < 0.01,'Head direction differs by >= 0.01 deg.');
assert(S.max_torso_axis_difference_deg < 0.01,'Torso direction differs by >= 0.01 deg.');
assert(S.max_head_origin_difference_m < 1e-4,'Head origin differs by >= 0.1 mm.');
assert(S.max_torso_origin_difference_m < 1e-4,'Torso origin differs by >= 0.1 mm.');
T=readtable(fullfile(workFolder,'poses.csv'));
assert(height(T)==S.frames && all(isfinite(T{:,2:end}),'all'));
fprintf('PASS: %s, %d frames; forward axes and origins agree with the production reference.\n', ...
    S.dataset,S.frames);
end
