function plotTriangulatedCheetah(P, varsUsed)
% P is [nBody x 3]
% varsUsed order:
% {'l_eye','r_eye','nose','l_ear','r_ear','spine','neck_base'}

    hold on

    idx.l_eye     = find(strcmp(varsUsed,'l_eye'));
    idx.r_eye     = find(strcmp(varsUsed,'r_eye'));
    idx.nose      = find(strcmp(varsUsed,'nose'));
    idx.l_ear     = find(strcmp(varsUsed,'l_ear'));
    idx.r_ear     = find(strcmp(varsUsed,'r_ear'));
    idx.spine     = find(strcmp(varsUsed,'spine'));
    idx.neck_base = find(strcmp(varsUsed,'neck_base'));

    % ----- points in the CRBR landmark colours -----
    if all(isfinite(P(idx.l_eye,:)))
        scatter3(P(idx.l_eye,1), P(idx.l_eye,2), P(idx.l_eye,3), 50, 'filled', 'MarkerFaceColor','b');
    end
    if all(isfinite(P(idx.r_eye,:)))
        scatter3(P(idx.r_eye,1), P(idx.r_eye,2), P(idx.r_eye,3), 50, 'filled', 'MarkerFaceColor','r');
    end
    if all(isfinite(P(idx.nose,:)))
        scatter3(P(idx.nose,1), P(idx.nose,2), P(idx.nose,3), 50, 'filled', 'MarkerFaceColor','k');
    end
    if all(isfinite(P(idx.neck_base,:)))
        scatter3(P(idx.neck_base,1), P(idx.neck_base,2), P(idx.neck_base,3), 50, 'filled', 'MarkerFaceColor','g');
    end
    if all(isfinite(P(idx.l_ear,:)))
        scatter3(P(idx.l_ear,1), P(idx.l_ear,2), P(idx.l_ear,3), 50, 'filled', 'MarkerFaceColor','c');
    end
    if all(isfinite(P(idx.r_ear,:)))
        scatter3(P(idx.r_ear,1), P(idx.r_ear,2), P(idx.r_ear,3), 50, 'filled', 'MarkerFaceColor','m');
    end
    if all(isfinite(P(idx.spine,:)))
        scatter3(P(idx.spine,1), P(idx.spine,2), P(idx.spine,3), 50, 'filled', 'MarkerFaceColor','y');
    end

    % ----- skeleton lines -----
    drawSeg(P, idx.l_eye, idx.r_eye);
    drawMidToPoint(P, idx.l_eye, idx.r_eye, idx.nose);
    drawSeg(P, idx.l_ear, idx.r_ear);
    drawSeg(P, idx.l_ear, idx.l_eye);
    drawSeg(P, idx.r_ear, idx.r_eye);
    drawSeg(P, idx.neck_base, idx.spine);
end

function drawSeg(P, ia, ib)
    if all(isfinite(P(ia,:))) && all(isfinite(P(ib,:)))
        plot3([P(ia,1) P(ib,1)], [P(ia,2) P(ib,2)], [P(ia,3) P(ib,3)], ...
            'k-', 'LineWidth', 2);
    end
end

function drawMidToPoint(P, ia, ib, ic)
    if all(isfinite(P(ia,:))) && all(isfinite(P(ib,:))) && all(isfinite(P(ic,:)))
        mid = (P(ia,:) + P(ib,:))/2;
        plot3([mid(1) P(ic,1)], [mid(2) P(ic,2)], [mid(3) P(ic,3)], ...
            'k-', 'LineWidth', 2);
    end
end
