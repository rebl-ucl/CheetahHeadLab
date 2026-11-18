function plotCheetah(pts)

hold on

l_eye     = pts.l_eye(:);
r_eye     = pts.r_eye(:);
nose      = pts.nose(:);
l_ear     = pts.l_ear(:);
r_ear     = pts.r_ear(:);
neck_base = pts.neck_base(:);
spine     = pts.spine(:);

scatter3(l_eye(1),l_eye(2),l_eye(3),50,'filled','MarkerFaceColor','b')
scatter3(r_eye(1),r_eye(2),r_eye(3),50,'filled','MarkerFaceColor','r')
scatter3(nose(1),nose(2),nose(3),50,'filled','MarkerFaceColor','k')
scatter3(l_ear(1),l_ear(2),l_ear(3),50,'filled','MarkerFaceColor','c')
scatter3(r_ear(1),r_ear(2),r_ear(3),50,'filled','MarkerFaceColor','m')
scatter3(neck_base(1),neck_base(2),neck_base(3),50,'filled','MarkerFaceColor','g')
scatter3(spine(1),spine(2),spine(3),50,'filled','MarkerFaceColor','y')

patch([l_eye(1), r_eye(1), nose(1)], ...
      [l_eye(2), r_eye(2), nose(2)], ...
      [l_eye(3), r_eye(3), nose(3)], ...
      'm', 'FaceAlpha',0.6)

% Ear triangle (l_ear, r_ear, nose) to visualise head orientation
patch([l_ear(1), r_ear(1), nose(1)], ...
      [l_ear(2), r_ear(2), nose(2)], ...
      [l_ear(3), r_ear(3), nose(3)], ...
      'c', 'FaceAlpha',0.4)

patch([l_eye(1), r_eye(1), neck_base(1)], ...
      [l_eye(2), r_eye(2), neck_base(2)], ...
      [l_eye(3), r_eye(3), neck_base(3)], ...
      'g', 'FaceAlpha',0.55)

plot3([neck_base(1), spine(1)], ...
      [neck_base(2), spine(2)], ...
      [neck_base(3), spine(3)], ...
      'k-', 'LineWidth', 2)

end
