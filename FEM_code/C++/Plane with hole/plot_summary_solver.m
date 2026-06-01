function plot_summary_solver(output_root)
% plot_summary_solver 绘制求解器耗时、迭代次数和 SOR omega 对比。

file = fullfile(output_root, 'summary_solver.csv');
if ~isfile(file)
    warning('没有找到 summary_solver.csv。');
    return;
end
T = readtable(file);
if isempty(T)
    warning('summary_solver.csv 为空。');
    return;
end

% 1. 直接法和迭代法总耗时对比
methods = string(T.method);
unique_methods = unique(methods, 'stable');
time_sum = zeros(numel(unique_methods), 1);
iter_mean = zeros(numel(unique_methods), 1);
for i = 1:numel(unique_methods)
    id = methods == unique_methods(i);
    time_sum(i) = sum(T.time_seconds(id));
    iter_mean(i) = mean(T.iterations(id));
end

fig = figure('Visible', 'off');
bar(time_sum);
grid on; box on;
set(gca, 'XTick', 1:numel(unique_methods), 'XTickLabel', unique_methods, 'XTickLabelRotation', 35);
ylabel('total time / s');
title('solver time compare');
saveas(fig, fullfile(output_root, 'solver_time_compare.png'));
close(fig);

fig = figure('Visible', 'off');
bar(iter_mean);
grid on; box on;
set(gca, 'XTick', 1:numel(unique_methods), 'XTickLabel', unique_methods, 'XTickLabelRotation', 35);
ylabel('mean iterations');
title('iterative solver iteration compare');
saveas(fig, fullfile(output_root, 'iterative_solver_iteration_compare.png'));
close(fig);

% 2. SOR omega 对比
is_sor = contains(methods, 'sor');
if any(is_sor)
    Ts = T(is_sor, :);
    fig = figure('Visible', 'off');
    scatter(Ts.omega, Ts.iterations, 50, Ts.time_seconds, 'filled');
    grid on; box on; colorbar;
    xlabel('omega'); ylabel('iterations');
    title('SOR omega compare: color = time / s');
    saveas(fig, fullfile(output_root, 'sor_omega_compare.png'));
    close(fig);
end

% 3. CG / PCG 对比
is_cg = contains(methods, 'cg');
if any(is_cg)
    Tc = T(is_cg, :);
    fig = figure('Visible', 'off');
    bar(categorical(string(Tc.method)), [Tc.iterations, Tc.time_seconds]);
    grid on; box on;
    ylabel('value');
    legend({'iterations', 'time / s'}, 'Location', 'best');
    title('CG and PCG compare');
    saveas(fig, fullfile(output_root, 'cg_pcg_compare.png'));
    close(fig);
end

% 4. dense / sparse 耗时对比
matrix_type = string(T.matrix_type);
fig = figure('Visible', 'off');
gscatter(1:height(T), T.time_seconds, matrix_type);
grid on; box on;
xlabel('solver record id'); ylabel('time / s');
title('dense sparse time compare');
saveas(fig, fullfile(output_root, 'dense_sparse_time_compare.png'));
close(fig);
end
