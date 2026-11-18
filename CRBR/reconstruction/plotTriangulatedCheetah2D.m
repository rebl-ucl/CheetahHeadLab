function plotTriangulatedCheetah2D(P, varsUsed, dims)
% varsUsed order:
% {'l_eye','r_eye','nose','l_ear','r_ear','spine','neck_base'}
    idx.l_eye     = find(strcmp(varsUsed,'l_eye'));
    idx.r_eye     = find(strcmp(varsUsed,'r_eye'));
    idx.nose      = find(strcmp(varsUsed,'nose'));
    idx.l_ear     = find(strcmp(varsUsed,'l_ear'));
    idx.r_ear     = find(strcmp(varsUsed,'r_ear'));
    idx.spine     = find(strcmp(varsUsed,'spine'));
    idx.neck_base = find(strcmp(varsUsed,'neck_base'));

    hold on

    scatterIfValid2D(P, idx.l_eye, dims, 'b');
    scatterIfValid2D(P, idx.r_eye, dims, 'r');
    scatterIfValid2D(P, idx.nose, dims, 'k');
    scatterIfValid2D(P, idx.neck_base, dims, 'g');
    scatterIfValid2D(P, idx.l_ear, dims, 'c');
    scatterIfValid2D(P, idx.r_ear, dims, 'm');
    scatterIfValid2D(P, idx.spine, dims, 'y');

    drawSeg2D(P, idx.l_eye, idx.r_eye, dims);
    drawMidToPoint2D(P, idx.l_eye, idx.r_eye, idx.nose, dims);
    drawSeg2D(P, idx.l_ear, idx.r_ear, dims);
    drawSeg2D(P, idx.l_ear, idx.l_eye, dims);
    drawSeg2D(P, idx.r_ear, idx.r_eye, dims);
    drawSeg2D(P, idx.neck_base, idx.spine, dims);
end

function scatterIfValid2D(P, i, dims, c)
    if all(isfinite(P(i,:)))
        scatter(P(i,dims(1)), P(i,dims(2)), 60, c, 'filled');
    end
end

function drawSeg2D(P, ia, ib, dims)
    if all(isfinite(P(ia,:))) && all(isfinite(P(ib,:)))
        plot([P(ia,dims(1)) P(ib,dims(1))], ...
             [P(ia,dims(2)) P(ib,dims(2))], ...
             'k-', 'LineWidth', 2);
    end
end

function drawMidToPoint2D(P, ia, ib, ic, dims)
    if all(isfinite(P(ia,:))) && all(isfinite(P(ib,:))) && all(isfinite(P(ic,:)))
        mid = (P(ia,:) + P(ib,:))/2;
        plot([mid(dims(1)) P(ic,dims(1))], ...
             [mid(dims(2)) P(ic,dims(2))], ...
             'k-', 'LineWidth', 2);
    end
end
