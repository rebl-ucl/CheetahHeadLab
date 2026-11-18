function [camUDx, camUDy] = undistortBlenderLabels(camLabels, intrin, nframes, varnames)

    camUDlabels = [];

    for i = 1:nframes
        rowlabels = [];

        for var = 1:numel(varnames)
            bp = string(varnames{var});
            xcol = bp + "_x";
            ycol = bp + "_y";

            if ismember(xcol, string(camLabels.Properties.VariableNames)) && ...
               ismember(ycol, string(camLabels.Properties.VariableNames))
                labelvalue = table2array(camLabels(i, [xcol, ycol]));  % 1x2
            else
                labelvalue = [NaN, NaN];
            end

            rowlabels = [rowlabels; labelvalue]; %#ok<AGROW>
        end

        % Undistort only finite points
        valid = all(isfinite(rowlabels), 2);
        udlabels = rowlabels;
        if any(valid)
            udlabels(valid,:) = undistortPoints(rowlabels(valid,:), intrin);
        end

        camUDlabels = [camUDlabels, udlabels]; %#ok<AGROW>
    end

    tpud = camUDlabels.';   % (2*nframes) x nParts

    % IMPORTANT: x/y are alternating ROWS here
    camUDx = array2table(tpud(1:2:end, :), 'VariableNames', varnames);
    camUDy = array2table(tpud(2:2:end, :), 'VariableNames', varnames);
end