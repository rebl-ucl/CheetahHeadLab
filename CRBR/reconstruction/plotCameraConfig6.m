
function plotCameraConfig6(cameraPosesArr, plotsize)
    cam1 = plotCamera(AbsolutePose=cameraPosesArr{1},Opacity=0,Color="blue",Size=plotsize, Label='1');
    grid on;axis equal;hold on;
    cam2 = plotCamera(AbsolutePose=cameraPosesArr{2},Opacity=0,Color="green",Size=plotsize, Label='2'); drawnow;
    cam3 = plotCamera(AbsolutePose=cameraPosesArr{3},Opacity=0,Color="magenta",Size=plotsize, Label='3'); drawnow;
    cam4 = plotCamera(AbsolutePose=cameraPosesArr{4},Opacity=0,Color="cyan",Size=plotsize, Label='4'); drawnow;
    cam5 = plotCamera(AbsolutePose=cameraPosesArr{5},Opacity=0,Color="red",Size=plotsize, Label='5'); drawnow;
    cam6 = plotCamera(AbsolutePose=cameraPosesArr{6},Opacity=0,Color="black",Size=plotsize, Label='6'); drawnow;
    % view([-0.205,-76.0743])
    title('6-Camera Concentric Setup')
    xlabel("X"); ylabel("Y"); zlabel("Z");
    hold off;
end