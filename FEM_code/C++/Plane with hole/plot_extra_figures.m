function plot_extra_figures(output_root)
% plot_extra_figures 生成大作业报告中建议补充的高级图像。
%
% 说明：
% 1. 本脚本只负责读取 C++ 输出的 csv/txt 文件并绘图，不参与有限元计算。
% 2. 若某些数据文件不存在，会自动跳过对应图像并给出 warning。
% 3. 推荐配合新版 C++ 文件使用，新版会额外输出 iterative_history.csv，
%    用于绘制残差-迭代步曲线和残差-时间曲线。

if nargin < 1 || isempty(output_root)
    output_root = fullfile('..', 'output');
end

if ~isfolder(output_root)
    warning('plot_extra_figures: 没有找到 output_root = %s', output_root);
    return;
end

ratio_list = find_ratio_list(output_root);
if isempty(ratio_list)
    warning('plot_extra_figures: 没有找到 ratio_* 工况文件夹。');
    return;
end

fprintf('开始绘制补充图像，共 %d 个 ratio 工况。\n', numel(ratio_list));

plot_residual_iteration(output_root, ratio_list);
plot_residual_time(output_root, ratio_list);
plot_solver_accuracy_time(output_root, ratio_list);
plot_sor_heatmap(output_root, ratio_list);
plot_best_omega(output_root, ratio_list);
plot_Kt_theory_error(output_root);
plot_solver_error_compare(output_root);
plot_tolerance_study(output_root);
plot_stiffness_spy(output_root, ratio_list);
plot_memory_compare(output_root, ratio_list);
plot_mesh_density_proxy(output_root, ratio_list);
plot_mesh_quality_box(output_root, ratio_list);
plot_hole_sigma_angle(output_root, ratio_list);
plot_hole_peak_zoom(output_root, ratio_list);
plot_curvature_Kt(output_root);
plot_total_load_check(output_root, ratio_list);
plot_constrained_dof_count(output_root, ratio_list);
plot_direct_pcg_difference(output_root, ratio_list);

fprintf('补充图像绘制完成。\n');
end

%% ============================================================
% 1. 残差-迭代步曲线
% ============================================================
function plot_residual_iteration(output_root, ratio_list)
fig = figure('Visible', 'off');
hold on; grid on; box on;
has_data = false;
for r = ratio_list(:)'
    hist_file = fullfile(output_root, sprintf('ratio_%g', r), 'iterative_history.csv');
    if ~isfile(hist_file), continue; end
    T = readtable(hist_file);
    if isempty(T), continue; end
    keys = unique(strcat(string(T.method), "_", string(T.matrix_type), "_", string(T.omega)), 'stable');
    for k = 1:numel(keys)
        id = strcat(string(T.method), "_", string(T.matrix_type), "_", string(T.omega)) == keys(k);
        Ti = T(id, :);
        if isempty(Ti), continue; end
        semilogy(Ti.iteration, max(Ti.residual, 1e-300), 'LineWidth', 1.1, ...
            'DisplayName', sprintf('b/a=%g %s', r, string(Ti.method(1))));
        has_data = true;
    end
end
xlabel('iteration');
ylabel('relative residual');
title('不同迭代法残差—迭代步曲线');
if has_data
    legend('Location', 'eastoutside', 'Interpreter', 'none');
    saveas(fig, fullfile(output_root, 'residual_iteration_curve.png'));
else
    warning('未找到 iterative_history.csv，跳过 residual_iteration_curve.png。');
end
close(fig);
end

%% ============================================================
% 2. 残差-时间曲线
% ============================================================
function plot_residual_time(output_root, ratio_list)
fig = figure('Visible', 'off');
hold on; grid on; box on;
has_data = false;
for r = ratio_list(:)'
    hist_file = fullfile(output_root, sprintf('ratio_%g', r), 'iterative_history.csv');
    if ~isfile(hist_file), continue; end
    T = readtable(hist_file);
    if isempty(T), continue; end
    keys = unique(strcat(string(T.method), "_", string(T.matrix_type), "_", string(T.omega)), 'stable');
    for k = 1:numel(keys)
        id = strcat(string(T.method), "_", string(T.matrix_type), "_", string(T.omega)) == keys(k);
        Ti = T(id, :);
        semilogy(Ti.time_seconds, max(Ti.residual, 1e-300), 'LineWidth', 1.1, ...
            'DisplayName', sprintf('b/a=%g %s', r, string(Ti.method(1))));
        has_data = true;
    end
end
xlabel('time / s');
ylabel('relative residual');
title('不同迭代法残差—时间曲线');
if has_data
    legend('Location', 'eastoutside', 'Interpreter', 'none');
    saveas(fig, fullfile(output_root, 'residual_time_curve.png'));
else
    warning('未找到 iterative_history.csv，跳过 residual_time_curve.png。');
end
close(fig);
end

%% ============================================================
% 3. 求解器精度—耗时散点图
% ============================================================
function plot_solver_accuracy_time(output_root, ratio_list)
S = collect_solver_summary(output_root, ratio_list);
if isempty(S), return; end
fig = figure('Visible', 'off');
scatter(S.time_seconds, max(S.final_residual, 1e-300), 70, S.ratio, 'filled');
set(gca, 'YScale', 'log');
grid on; box on; colorbar;
xlabel('time / s');
ylabel('final relative residual');
title('求解器精度—耗时散点图：颜色表示 b/a');
text(S.time_seconds, max(S.final_residual, 1e-300), string(S.method), ...
    'FontSize', 7, 'Interpreter', 'none');
saveas(fig, fullfile(output_root, 'solver_accuracy_time_scatter.png'));
close(fig);
end

%% ============================================================
% 4. SOR omega-ratio-iterations 热力图
% ============================================================
function plot_sor_heatmap(output_root, ratio_list)
S = collect_solver_summary(output_root, ratio_list);
if isempty(S), return; end
id = contains(string(S.method), 'sor');
if ~any(id)
    warning('没有 SOR 求解器记录，跳过 SOR 热力图。若需要该图，请关闭 auto_skip_expensive_solver 或用较小网格运行 SOR。');
    return;
end
Ts = S(id, :);
omegas = unique(Ts.omega);
ratios = unique(Ts.ratio);
Z = nan(numel(ratios), numel(omegas));
for i = 1:numel(ratios)
    for j = 1:numel(omegas)
        k = Ts.ratio == ratios(i) & abs(Ts.omega - omegas(j)) < 1e-12;
        if any(k)
            Z(i, j) = min(Ts.iterations(k));
        end
    end
end
fig = figure('Visible', 'off');
imagesc(omegas, ratios, Z);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('SOR \omega'); ylabel('b/a');
title('SOR omega-ratio-iterations 热力图');
saveas(fig, fullfile(output_root, 'sor_omega_ratio_iterations_heatmap.png'));
close(fig);
end

%% ============================================================
% 5. 最优 omega 随 b/a 变化图
% ============================================================
function plot_best_omega(output_root, ratio_list)
S = collect_solver_summary(output_root, ratio_list);
if isempty(S), return; end
id = contains(string(S.method), 'sor') & S.converged == 1;
if ~any(id)
    warning('没有收敛 SOR 记录，跳过 best_omega_vs_ratio.png。');
    return;
end
Ts = S(id, :);
ratios = unique(Ts.ratio);
best_omega = nan(size(ratios));
best_iter = nan(size(ratios));
for i = 1:numel(ratios)
    k = find(Ts.ratio == ratios(i));
    [best_iter(i), loc] = min(Ts.iterations(k));
    best_omega(i) = Ts.omega(k(loc));
end
fig = figure('Visible', 'off');
yyaxis left;
plot(ratios, best_omega, '-o', 'LineWidth', 1.4);
ylabel('best \omega');
yyaxis right;
plot(ratios, best_iter, '-s', 'LineWidth', 1.4);
ylabel('iterations at best \omega');
grid on; box on;
xlabel('b/a');
title('最优 SOR \omega 随 b/a 变化');
saveas(fig, fullfile(output_root, 'best_omega_vs_ratio.png'));
close(fig);
end

%% ============================================================
% 6. Kt_FE、Kt_theory、error_Kt 三合一图
% ============================================================
function plot_Kt_theory_error(output_root)
file = fullfile(output_root, 'summary_ratio.csv');
if ~isfile(file), return; end
T = readtable(file);
if isempty(T), return; end
Kt_FE = T.Kt_x;
Kt_theory = 1 + 2 * T.ratio;  % Inglis 椭圆孔无限板参考解，有限板只作趋势参照。
error_Kt = abs(Kt_FE - Kt_theory) ./ max(abs(Kt_theory), 1e-30);
fig = figure('Visible', 'off');
yyaxis left;
plot(T.ratio, Kt_FE, '-o', 'LineWidth', 1.5); hold on;
plot(T.ratio, Kt_theory, '--s', 'LineWidth', 1.5);
ylabel('K_t');
yyaxis right;
plot(T.ratio, error_Kt, '-.^', 'LineWidth', 1.5);
ylabel('relative error');
grid on; box on;
xlabel('b/a');
title('Kt\_FE、Kt\_theory、error\_Kt 三合一图');
legend({'Kt\_FE by \sigma_x', 'Kt\_theory = 1 + 2b/a', 'relative error'}, 'Location', 'best');
saveas(fig, fullfile(output_root, 'Kt_FE_theory_error.png'));
close(fig);

Tout = table(T.ratio, Kt_FE, Kt_theory, error_Kt, ...
    'VariableNames', {'ratio','Kt_FE_sigma_x','Kt_theory_Inglis_reference','relative_error'});
writetable(Tout, fullfile(output_root, 'Kt_theory_error.csv'));
end

%% ============================================================
% 7. 不同求解器 error 对比图
% ============================================================
function plot_solver_error_compare(output_root)
file = fullfile(output_root, 'compare_solver_result.csv');
if ~isfile(file), return; end
T = readtable(file);
if isempty(T), return; end
fig = figure('Visible', 'off');
scatter(T.ratio, max(T.error_residual, 1e-300), 60, categorical(string(T.method_compare)), 'filled');
set(gca, 'YScale', 'log');
grid on; box on;
xlabel('b/a');
ylabel('relative residual difference');
title('不同求解器误差对比图');
legend('Location', 'eastoutside', 'Interpreter', 'none');
saveas(fig, fullfile(output_root, 'solver_error_compare.png'));
close(fig);
end

%% ============================================================
% 8. tol 对 error_Kt 和耗时影响图
% ============================================================
function plot_tolerance_study(output_root)
file = fullfile(output_root, 'tolerance_study.csv');
if ~isfile(file)
    warning('没有 tolerance_study.csv，跳过 tol_error_Kt_time.png。需要多次修改 tol 运行后手动汇总。');
    return;
end
T = readtable(file);
need = {'tol','error_Kt','time_seconds'};
if ~all(ismember(need, T.Properties.VariableNames))
    warning('tolerance_study.csv 缺少 tol/error_Kt/time_seconds 字段。');
    return;
end
fig = figure('Visible', 'off');
yyaxis left;
loglog(T.tol, T.error_Kt, '-o', 'LineWidth', 1.4);
ylabel('error\_Kt');
yyaxis right;
semilogx(T.tol, T.time_seconds, '-s', 'LineWidth', 1.4);
ylabel('time / s');
grid on; box on;
xlabel('tolerance');
title('tol 对 error\_Kt 和耗时的影响');
saveas(fig, fullfile(output_root, 'tol_error_Kt_time.png'));
close(fig);
end

%% ============================================================
% 9. 总体刚度矩阵稀疏结构 spy 图
% ============================================================
function plot_stiffness_spy(output_root, ratio_list)
% 优先选择 ratio=3，否则选第一个工况。
r = 3;
if ~ismember(r, ratio_list), r = ratio_list(1); end
case_dir = fullfile(output_root, sprintf('ratio_%g', r));
node_file = fullfile(case_dir, 'nodes.csv');
elem_file = fullfile(case_dir, 'elements.csv');
if ~isfile(node_file) || ~isfile(elem_file), return; end
nodes = readtable(node_file);
elems = readtable(elem_file);
n_dof = 2 * height(nodes);
I = [];
J = [];
E = [elems.n1, elems.n2, elems.n3, elems.n4] + 1;
for e = 1:size(E,1)
    dof = zeros(1,8);
    for k = 1:4
        nid = E(e,k) - 1;
        dof(2*k-1) = 2*nid + 1;
        dof(2*k)   = 2*nid + 2;
    end
    [ii,jj] = ndgrid(dof,dof);
    I = [I; ii(:)]; %#ok<AGROW>
    J = [J; jj(:)]; %#ok<AGROW>
end
Kpat = sparse(I, J, true, n_dof, n_dof);
fig = figure('Visible', 'off');
spy(Kpat, 2);
title(sprintf('总体刚度矩阵 sparse 结构图，b/a=%g', r));
xlabel('DOF index'); ylabel('DOF index');
saveas(fig, fullfile(output_root, sprintf('stiffness_spy_ratio_%g.png', r)));
close(fig);
end

%% ============================================================
% 10. dense/sparse 内存占用对比图
% ============================================================
function plot_memory_compare(output_root, ratio_list)
ratio = [];
dense_MB = [];
sparse_MB = [];
for r = ratio_list(:)'
    info = read_check_info(fullfile(output_root, sprintf('ratio_%g', r), 'check_info.txt'));
    if isempty(info), continue; end
    if ~isfield(info, 'n_dof') || ~isfield(info, 'sparse_nnz'), continue; end
    n = info.n_dof;
    nnzK = info.sparse_nnz;
    dense_bytes = double(n) * double(n) * 8;
    csr_bytes = double(nnzK) * (8 + 4) + double(n + 1) * 4;
    ratio(end+1,1) = r; %#ok<AGROW>
    dense_MB(end+1,1) = dense_bytes / 1024^2; %#ok<AGROW>
    sparse_MB(end+1,1) = csr_bytes / 1024^2; %#ok<AGROW>
end
if isempty(ratio), return; end
fig = figure('Visible', 'off');
bar(ratio, [dense_MB, sparse_MB]);
grid on; box on;
xlabel('b/a'); ylabel('estimated memory / MB');
legend({'dense K', 'CSR sparse K'}, 'Location', 'best');
title('dense/sparse 内存占用对比图');
saveas(fig, fullfile(output_root, 'dense_sparse_memory_compare.png'));
close(fig);
end

%% ============================================================
% 11. 网格规模—error_Kt—time 图
% ============================================================
function plot_mesh_density_proxy(output_root, ratio_list)
Tsum_file = fullfile(output_root, 'summary_ratio.csv');
if ~isfile(Tsum_file), return; end
Tsum = readtable(Tsum_file);
Tkt = readtable(fullfile(output_root, 'Kt_theory_error.csv'));
ratio = [];
n_elem = [];
time_total = [];
error_Kt = [];
for r = ratio_list(:)'
    info = read_check_info(fullfile(output_root, sprintf('ratio_%g', r), 'check_info.txt'));
    time_file = fullfile(output_root, sprintf('ratio_%g', r), 'time_summary.csv');
    if isempty(info) || ~isfile(time_file), continue; end
    Ttime = readtable(time_file);
    idt = strcmp(string(Ttime.step_name), 'total_case_time');
    ratio(end+1,1) = r; %#ok<AGROW>
    n_elem(end+1,1) = info.n_elem; %#ok<AGROW>
    if any(idt), time_total(end+1,1) = Ttime.time_seconds(find(idt,1)); else, time_total(end+1,1) = nan; end %#ok<AGROW>
    idk = abs(Tkt.ratio - r) < 1e-12;
    if any(idk), error_Kt(end+1,1) = Tkt.relative_error(find(idk,1)); else, error_Kt(end+1,1) = nan; end %#ok<AGROW>
end
if isempty(ratio), return; end
fig = figure('Visible', 'off');
scatter3(n_elem, time_total, error_Kt, 80, ratio, 'filled');
grid on; box on; colorbar;
xlabel('n\_elem'); ylabel('total time / s'); zlabel('error\_Kt');
title('网格规模—error\_Kt—time 图');
saveas(fig, fullfile(output_root, 'mesh_density_error_Kt_time.png'));
close(fig);
end

%% ============================================================
% 12. detJ 或 aspect ratio 箱线图
% ============================================================
function plot_mesh_quality_box(output_root, ratio_list)
all_ratio = [];
all_ar = [];
all_detJ = [];
for r = ratio_list(:)'
    case_dir = fullfile(output_root, sprintf('ratio_%g', r));
    node_file = fullfile(case_dir, 'nodes.csv');
    elem_file = fullfile(case_dir, 'elements.csv');
    if ~isfile(node_file) || ~isfile(elem_file), continue; end
    nodes = readtable(node_file);
    elems = readtable(elem_file);
    [ar, detJ] = calculate_mesh_quality_arrays(nodes, elems);
    all_ratio = [all_ratio; r * ones(numel(ar),1)]; %#ok<AGROW>
    all_ar = [all_ar; ar(:)]; %#ok<AGROW>
    all_detJ = [all_detJ; detJ(:)]; %#ok<AGROW>
end
if isempty(all_ar), return; end
fig = figure('Visible', 'off');
boxplot(all_ar, all_ratio);
grid on; box on;
xlabel('b/a'); ylabel('aspect ratio');
title('aspect ratio 箱线图');
saveas(fig, fullfile(output_root, 'aspect_ratio_boxplot.png'));
close(fig);
fig = figure('Visible', 'off');
boxplot(all_detJ, all_ratio);
grid on; box on;
xlabel('b/a'); ylabel('detJ at Gauss points');
title('detJ 箱线图');
saveas(fig, fullfile(output_root, 'detJ_boxplot.png'));
close(fig);
end

%% ============================================================
% 13. 孔边 sigma_x 沿角度分布曲线
% ============================================================
function plot_hole_sigma_angle(output_root, ratio_list)
fig_all = figure('Visible', 'off');
hold on; grid on; box on;
has_any = false;
for r = ratio_list(:)'
    case_dir = fullfile(output_root, sprintf('ratio_%g', r));
    node_file = fullfile(case_dir, 'nodes.csv');
    elem_file = fullfile(case_dir, 'elements.csv');
    stress_file = fullfile(case_dir, 'stress.csv');
    sum_file = fullfile(case_dir, 'case_summary.csv');
    if ~all(cellfun(@isfile, {node_file, elem_file, stress_file, sum_file})), continue; end
    nodes = readtable(node_file);
    elems = readtable(elem_file);
    stress = readtable(stress_file);
    S = readtable(sum_file);
    E = [elems.n1, elems.n2, elems.n3, elems.n4] + 1;
    hole_flag = nodes.flag_hole;
    is_hole_elem = any(hole_flag(E) == 1, 2);
    idx = find(is_hole_elem);
    if isempty(idx), continue; end
    x = stress.x_center(idx);
    y = stress.y_center(idx);
    theta = atan2(y, x - 42.0);
    theta(theta < 0) = theta(theta < 0) + pi;
    [theta, order] = sort(theta);
    sx = stress.sigma_x(idx(order));
    fig = figure('Visible', 'off');
    plot(theta * 180/pi, sx, '-o', 'LineWidth', 1.2, 'MarkerSize', 3);
    grid on; box on;
    xlabel('\theta / degree'); ylabel('\sigma_x / MPa');
    title(sprintf('孔边 \sigma_x 沿角度分布曲线，b/a=%g', r));
    saveas(fig, fullfile(case_dir, 'hole_sigma_x_angle.png'));
    close(fig);
    plot(theta * 180/pi, sx, 'LineWidth', 1.2, 'DisplayName', sprintf('b/a=%g', r));
    has_any = true;
end
if has_any
    xlabel('\theta / degree'); ylabel('\sigma_x / MPa');
    title('孔边 \sigma_x 沿角度分布曲线对比');
    legend('Location', 'best');
    saveas(fig_all, fullfile(output_root, 'hole_sigma_x_angle_compare.png'));
else
    warning('未找到孔边单元，跳过孔边 sigma_x 曲线。');
end
close(fig_all);
end

%% ============================================================
% 14. 孔边峰值区域局部放大图
% ============================================================
function plot_hole_peak_zoom(output_root, ratio_list)
for r = ratio_list(:)'
    case_dir = fullfile(output_root, sprintf('ratio_%g', r));
    node_file = fullfile(case_dir, 'nodes.csv');
    elem_file = fullfile(case_dir, 'elements.csv');
    stress_file = fullfile(case_dir, 'stress.csv');
    if ~all(cellfun(@isfile, {node_file, elem_file, stress_file})), continue; end
    nodes = readtable(node_file);
    elems = readtable(elem_file);
    stress = readtable(stress_file);
    [~, id] = max(stress.sigma_vm);
    xc = stress.x_center(id); yc = stress.y_center(id);
    V = [nodes.x, nodes.y];
    F = [elems.n1, elems.n2, elems.n3, elems.n4] + 1;
    fig = figure('Visible', 'off');
    patch('Faces', F, 'Vertices', V, 'FaceVertexCData', stress.sigma_vm, ...
        'FaceColor', 'flat', 'EdgeColor', [0.6 0.6 0.6]);
    hold on;
    plot(xc, yc, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 12);
    axis equal; box on; colorbar; colormap(jet);
    xlim([xc - 8, xc + 8]); ylim([max(0, yc - 8), yc + 8]);
    xlabel('x / mm'); ylabel('y / mm');
    title(sprintf('孔边峰值区域局部放大，b/a=%g', r));
    saveas(fig, fullfile(case_dir, 'hole_peak_zoom.png'));
    close(fig);
end
end

%% ============================================================
% 15. 曲率半径 rho 与 Kt_FE 关系图
% ============================================================
function plot_curvature_Kt(output_root)
file = fullfile(output_root, 'summary_ratio.csv');
if ~isfile(file), return; end
T = readtable(file);
if isempty(T), return; end
rho = zeros(height(T),1);
for i = 1:height(T)
    theta = atan2(T.max_sigma_vm_y(i), T.max_sigma_vm_x(i) - 42.0);
    a = T.a(i); b = T.b(i);
    rho(i) = ((a^2 * sin(theta)^2 + b^2 * cos(theta)^2)^(3/2)) / max(a*b, 1e-30);
end
fig = figure('Visible', 'off');
scatter(rho, T.Kt_x, 80, T.ratio, 'filled');
grid on; box on; colorbar;
xlabel('local curvature radius \rho / mm');
ylabel('Kt\_FE by \sigma_x');
title('曲率半径 \rho 与 Kt\_FE 关系图');
saveas(fig, fullfile(output_root, 'curvature_radius_vs_Kt_FE.png'));
close(fig);
end

%% ============================================================
% 16. 实际总载荷与理论总载荷对比图
% ============================================================
function plot_total_load_check(output_root, ratio_list)
ratio = [];
actual = [];
theory = [];
for r = ratio_list(:)'
    info = read_check_info(fullfile(output_root, sprintf('ratio_%g', r), 'check_info.txt'));
    if isempty(info), continue; end
    ratio(end+1,1) = r; %#ok<AGROW>
    actual(end+1,1) = getfield_safe(info, 'actual_total_force', nan); %#ok<GFLD,AGROW>
    theory(end+1,1) = getfield_safe(info, 'theory_total_force', nan); %#ok<GFLD,AGROW>
end
if isempty(ratio), return; end
fig = figure('Visible', 'off');
plot(ratio, actual, '-o', 'LineWidth', 1.4); hold on;
plot(ratio, theory, '--s', 'LineWidth', 1.4);
grid on; box on;
xlabel('b/a'); ylabel('total force / N');
title('实际总载荷与理论总载荷对比图');
legend({'actual total load', 'theory total load'}, 'Location', 'best');
saveas(fig, fullfile(output_root, 'total_load_actual_theory_compare.png'));
close(fig);
end

%% ============================================================
% 17. 约束自由度数量统计图
% ============================================================
function plot_constrained_dof_count(output_root, ratio_list)
ratio = [];
left = [];
sym = [];
for r = ratio_list(:)'
    info = read_check_info(fullfile(output_root, sprintf('ratio_%g', r), 'check_info.txt'));
    if isempty(info), continue; end
    ratio(end+1,1) = r; %#ok<AGROW>
    left(end+1,1) = getfield_safe(info, 'left_fixed_dof', nan); %#ok<GFLD,AGROW>
    sym(end+1,1) = getfield_safe(info, 'symmetry_fixed_dof', nan); %#ok<GFLD,AGROW>
end
if isempty(ratio), return; end
fig = figure('Visible', 'off');
bar(ratio, [left, sym], 'stacked');
grid on; box on;
xlabel('b/a'); ylabel('constrained DOF count');
legend({'left fixed DOF', 'symmetry fixed DOF'}, 'Location', 'best');
title('约束自由度数量统计图');
saveas(fig, fullfile(output_root, 'constrained_dof_count.png'));
close(fig);
end

%% ============================================================
% 18. 直接法与 PCG 位移差异云图
% ============================================================
function plot_direct_pcg_difference(output_root, ratio_list)
% 需要额外输出 displacement_direct.csv 和 displacement_pcg.csv。
% 如果没有这两个文件，本脚本会跳过该图。
for r = ratio_list(:)'
    case_dir = fullfile(output_root, sprintf('ratio_%g', r));
    f1 = fullfile(case_dir, 'displacement_direct.csv');
    f2 = fullfile(case_dir, 'displacement_pcg.csv');
    elem_file = fullfile(case_dir, 'elements.csv');
    if ~all(cellfun(@isfile, {f1, f2, elem_file})), continue; end
    D1 = readtable(f1); D2 = readtable(f2); elems = readtable(elem_file);
    du = D2.ux - D1.ux;
    dv = D2.uy - D1.uy;
    err = sqrt(du.^2 + dv.^2);
    V = [D1.x, D1.y];
    F = [elems.n1, elems.n2, elems.n3, elems.n4] + 1;
    fig = figure('Visible', 'off');
    patch('Faces', F, 'Vertices', V, 'FaceVertexCData', err, ...
        'FaceColor', 'interp', 'EdgeColor', 'none');
    axis equal; box on; colorbar; colormap(jet);
    xlabel('x / mm'); ylabel('y / mm');
    title(sprintf('直接法与 PCG 位移差异云图，b/a=%g', r));
    saveas(fig, fullfile(case_dir, 'direct_pcg_displacement_difference.png'));
    close(fig);
end
end

%% ============================================================
% 公共辅助函数
% ============================================================
function ratio_list = find_ratio_list(output_root)
d = dir(fullfile(output_root, 'ratio_*'));
ratio_list = [];
for i = 1:numel(d)
    if ~d(i).isdir, continue; end
    token = regexp(d(i).name, '^ratio_([0-9\.]+)$', 'tokens');
    if isempty(token), continue; end
    ratio_list(end+1) = str2double(token{1}{1}); %#ok<AGROW>
end
ratio_list = sort(ratio_list);
end

function S = collect_solver_summary(output_root, ratio_list)
S = table();
for r = ratio_list(:)'
    f = fullfile(output_root, sprintf('ratio_%g', r), 'iterative_solver_summary.csv');
    if ~isfile(f), continue; end
    T = readtable(f);
    if isempty(T), continue; end
    T.ratio = r * ones(height(T),1);
    S = [S; T]; %#ok<AGROW>
    fd = fullfile(output_root, sprintf('ratio_%g', r), 'direct_solver_summary.csv');
    if isfile(fd)
        D = readtable(fd);
        if ~isempty(D)
            D.omega = zeros(height(D),1);
            D.converged = ones(height(D),1);
            D.iterations = ones(height(D),1);
            D.ratio = r * ones(height(D),1);
            D = D(:, T.Properties.VariableNames);
            S = [S; D]; %#ok<AGROW>
        end
    end
end
if isempty(S)
    warning('没有找到求解器 summary 文件。');
end
end

function info = read_check_info(file)
info = struct();
if ~isfile(file)
    info = [];
    return;
end
fid = fopen(file, 'r');
if fid < 0
    info = [];
    return;
end
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    parts = regexp(line, '^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([^\n]+)\s*$', 'tokens');
    if isempty(parts), continue; end
    key = parts{1}{1};
    valstr = strtrim(parts{1}{2});
    val = str2double(regexp(valstr, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match', 'once'));
    if ~isnan(val)
        info.(key) = val;
    end
end
fclose(fid);
end

function v = getfield_safe(s, name, default_v)
if isfield(s, name), v = s.(name); else, v = default_v; end
end

function [ar, detJ_all] = calculate_mesh_quality_arrays(nodes, elems)
V = [nodes.x, nodes.y];
E = [elems.n1, elems.n2, elems.n3, elems.n4] + 1;
n_elem = size(E,1);
ar = zeros(n_elem,1);
detJ_all = [];
g = 1 / sqrt(3);
gps = [-g,-g; g,-g; g,g; -g,g];
for e = 1:n_elem
    xy = V(E(e,:), :);
    len = [norm(xy(2,:)-xy(1,:)), norm(xy(3,:)-xy(2,:)), ...
           norm(xy(4,:)-xy(3,:)), norm(xy(1,:)-xy(4,:))];
    ar(e) = max(len) / max(min(len), 1e-30);
    for k = 1:4
        xi = gps(k,1); eta = gps(k,2);
        dN_dxi = 0.25 * [-(1-eta); (1-eta); (1+eta); -(1+eta)];
        dN_deta = 0.25 * [-(1-xi); -(1+xi); (1+xi); (1-xi)];
        J = [dN_dxi' * xy(:,1), dN_dxi' * xy(:,2); ...
             dN_deta' * xy(:,1), dN_deta' * xy(:,2)];
        detJ_all(end+1,1) = det(J); %#ok<AGROW>
    end
end
end
