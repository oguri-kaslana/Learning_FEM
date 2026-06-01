function plot_summary_ratio(output_root)
% plot_summary_ratio 绘制不同 b/a 下最大应力和应力集中系数变化。

file = fullfile(output_root, 'summary_ratio.csv');
if ~isfile(file)
    warning('没有找到 summary_ratio.csv。');
    return;
end
T = readtable(file);

fig = figure('Visible', 'off');
plot(T.ratio, T.max_sigma_vm, '-o', 'LineWidth', 1.5);
grid on; box on;
xlabel('b/a'); ylabel('max von Mises stress / MPa');
title('max\_sigma\_vm vs b/a');
saveas(fig, fullfile(output_root, 'max_sigma_vm_vs_ratio.png'));
close(fig);

fig = figure('Visible', 'off');
plot(T.ratio, T.Kt_vm, '-o', 'LineWidth', 1.5);
grid on; box on;
xlabel('b/a'); ylabel('K_t by von Mises');
title('Kt\_vm vs b/a');
saveas(fig, fullfile(output_root, 'Kt_vm_vs_ratio.png'));
close(fig);

fig = figure('Visible', 'off');
plot(T.ratio, T.max_umag, '-o', 'LineWidth', 1.5);
grid on; box on;
xlabel('b/a'); ylabel('max displacement magnitude / mm');
title('max displacement vs b/a');
saveas(fig, fullfile(output_root, 'max_disp_vs_ratio.png'));
close(fig);
end
