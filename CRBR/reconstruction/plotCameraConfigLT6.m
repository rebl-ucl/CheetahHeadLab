function plotCameraConfigLT6(cameraPosesArr, plotsize, camIDs)
    ncams = numel(camIDs);

    plotCamera(AbsolutePose=cameraPosesArr{1}, Opacity=0, Color="blue", ...
        Size=plotsize, Label=num2str(camIDs(1)));
    grid on; axis equal; hold on;

    plotCamera(AbsolutePose=cameraPosesArr{2}, Opacity=0, Color="green", ...
        Size=plotsize, Label=num2str(camIDs(2)));
    drawnow;

    switch ncams
        case 3
            plotCamera(AbsolutePose=cameraPosesArr{3}, Opacity=0, Color="magenta", ...
                Size=plotsize, Label=num2str(camIDs(3)));
            drawnow;

        case 4
            plotCamera(AbsolutePose=cameraPosesArr{3}, Opacity=0, Color="magenta", ...
                Size=plotsize, Label=num2str(camIDs(3)));
            drawnow;
            plotCamera(AbsolutePose=cameraPosesArr{4}, Opacity=0, Color="cyan", ...
                Size=plotsize, Label=num2str(camIDs(4)));
            drawnow;

        case 5
            plotCamera(AbsolutePose=cameraPosesArr{3}, Opacity=0, Color="magenta", ...
                Size=plotsize, Label=num2str(camIDs(3)));
            drawnow;
            plotCamera(AbsolutePose=cameraPosesArr{4}, Opacity=0, Color="cyan", ...
                Size=plotsize, Label=num2str(camIDs(4)));
            drawnow;
            plotCamera(AbsolutePose=cameraPosesArr{5}, Opacity=0, Color="red", ...
                Size=plotsize, Label=num2str(camIDs(5)));
            drawnow;
    end

    title(strcat(num2str(ncams), ' Camera Concentric Setup'))
    xlabel("X"); ylabel("Y"); zlabel("Z");
    hold off;
end