function plot_mesh(case_dir)
% plot_mesh 绘制 Q4 网格图。

nodes_file = fullfile(case_dir, 'nodes.csv');
elems_file = fullfile(case_dir, 'elements.csv');
summary_file = fullfile(case_dir, 'case_summary.csv');

nodes = readtable(nodes_file);
elems = readtable(elems_file);
summary = readtable(summary_file);

V = [nodes.x, nodes.y];
F = [elems.n1, elems.n2, elems.n3, elems.n4] + 1;

fig = figure('Visible', 'off');
patch('Faces', F, 'Vertices', V, ...
      'FaceColor', 'none', 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.5);
axis equal; grid on; box on;
xlabel('x / mm'); ylabel('y / mm');
title(sprintf('Q4 mesh: b/a = %.2f, a = %.3f mm, b = %.3f mm', ...
      summary.ratio(1), summary.a(1), summary.b(1)));

saveas(fig, fullfile(case_dir, 'mesh.png'));
close(fig);
end
