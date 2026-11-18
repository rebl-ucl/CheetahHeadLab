function r = ASseedResidualFromModel( ...
    p, seedAS, Rz, Rx, Ry, ...
    scaleHead, scaleTorsoLength)

    [~, ~, mdl] = cheetah3D( ...
        p, Rz, Rx, Ry, ...
        scaleHead, scaleTorsoLength);

    % Ears are only seeded when the seed source (e.g. multi-view
    % triangulation) provides them; the AS seed has no ear positions, so
    % skip any landmark missing from seedAS.
    names = {'l_eye','r_eye','nose','l_ear','r_ear','spine','neck_base'};

    r = [];
    for k = 1:numel(names)
        nm = names{k};
        if ~isfield(seedAS, nm)
            continue;
        end
        p_model = mdl.(nm)(:);
        p_seed  = seedAS.(nm)(:);
        r = [r; p_model - p_seed];
    end
end
