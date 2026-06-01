function plot_stress(case_dir, stress_name)
% plot_stress 绘制单元中心应力云图。
% stress_name 可取：sigma_x, sigma_y, tau_xy, sigma_vm。

nodes = readtable(fullfile(case_dir, 'nodes.csv'));
elems = readtable(fullfile(case_dir, 'elements.csv'));
stress = readtable(fullfile(case_dir, 'stress.csv'));
summary = readtable(fullfile(case_dir, 'case_summary.csv'));

V = [nodes.x, nodes.y];
F = [elems.n1, elems.n2, elems.n3, elems.n4] + 1;

if ~ismember(stress_name, stress.Properties.VariableNames)
    error('stress.csv 中不存在字段 %s。', stress_name);
end

C = stress.(stress_name);

fig = figure('Visible', 'off');
patch('Faces', F, 'Vertices', V, ...
      'FaceVertexCData', C, 'FaceColor', 'flat', ...
      'EdgeColor', 'none');
axis equal; box on;
xlabel('x / mm'); ylabel('y / mm');
colorbar;
title(sprintf('%s: b/a = %.2f, a = %.3f mm, b = %.3f mm', ...
      strrep(stress_name, '_', '\_'), summary.ratio(1), summary.a(1), summary.b(1)));

colormap(jet);
saveas(fig, fullfile(case_dir, [stress_name, '.png']));
close(fig);
end
