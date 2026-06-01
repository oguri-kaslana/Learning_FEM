function plot_deformation(case_dir, scale)
% plot_deformation 绘制变形图。
% scale 为位移放大系数。

if nargin < 2
    scale = 20;
end

nodes = readtable(fullfile(case_dir, 'nodes.csv'));
elems = readtable(fullfile(case_dir, 'elements.csv'));
disp_data = readtable(fullfile(case_dir, 'displacement.csv'));
summary = readtable(fullfile(case_dir, 'case_summary.csv'));

V0 = [nodes.x, nodes.y];
V1 = [nodes.x + scale * disp_data.ux, nodes.y + scale * disp_data.uy];
F = [elems.n1, elems.n2, elems.n3, elems.n4] + 1;

fig = figure('Visible', 'off');
patch('Faces', F, 'Vertices', V0, ...
      'FaceColor', 'none', 'EdgeColor', [0.75 0.75 0.75], 'LineWidth', 0.5);
hold on;
patch('Faces', F, 'Vertices', V1, ...
      'FaceColor', 'none', 'EdgeColor', [0 0.2 0.8], 'LineWidth', 0.7);
axis equal; grid on; box on;
xlabel('x / mm'); ylabel('y / mm');
title(sprintf('Deformation: b/a = %.2f, scale = %.1f', summary.ratio(1), scale));
legend({'undeformed', 'deformed'}, 'Location', 'best');

saveas(fig, fullfile(case_dir, 'deformation.png'));
close(fig);
end
