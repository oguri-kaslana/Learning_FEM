function plot_summary_time(output_root)
% plot_summary_time 绘制步骤耗时条形图，并绘制基础版/加速版加速比。

file = fullfile(output_root, 'summary_time.csv');
if isfile(file)
    T = readtable(file);
    if ~isempty(T)
        [time_sorted, id] = sort(T.time_seconds, 'descend');
        n_show = min(25, numel(id));
        fig = figure('Visible', 'off');
        bar(time_sorted(1:n_show));
        grid on; box on;
        set(gca, 'XTick', 1:n_show, 'XTickLabel', string(T.step_name(id(1:n_show))), 'XTickLabelRotation', 60);
        ylabel('time / s');
        title('top step time cost');
        saveas(fig, fullfile(output_root, 'step_time_bar.png'));
        close(fig);
    end
else
    warning('没有找到 summary_time.csv。');
end

compare_file = fullfile(output_root, 'compare_base_fast.csv');
if isfile(compare_file)
    C = readtable(compare_file);
    if ~isempty(C)
        fig = figure('Visible', 'off');
        plot(C.ratio, C.speedup, '-o', 'LineWidth', 1.5);
        grid on; box on;
        xlabel('b/a'); ylabel('speedup = time_base / time_fast');
        title('base fast speedup compare');
        saveas(fig, fullfile(output_root, 'base_fast_speedup_compare.png'));
        close(fig);
    end
end
end
