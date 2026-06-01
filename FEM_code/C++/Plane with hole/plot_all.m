%% plot_all.m
% 读取 C++ 输出结果并批量绘制所有云图和汇总图。
% MATLAB 只负责可视化，不参与有限元计算。

clear; clc; close all;

output_root = fullfile('..', 'output');
ratio_list = [1, 2, 3, 4, 5];
deformation_scale = 20;

if ~isfolder(output_root)
    error('没有找到 output 文件夹。请先运行 C++ 程序生成结果。');
end

for i = 1:numel(ratio_list)
    ratio = ratio_list(i);
    case_dir = fullfile(output_root, sprintf('ratio_%g', ratio));
    if ~isfolder(case_dir)
        warning('跳过 %s：文件夹不存在。', case_dir);
        continue;
    end

    fprintf('Plot ratio = %g\n', ratio);
    plot_mesh(case_dir);
    plot_deformation(case_dir, deformation_scale);
    plot_stress(case_dir, 'sigma_x');
    plot_stress(case_dir, 'sigma_y');
    plot_stress(case_dir, 'tau_xy');
    plot_stress(case_dir, 'sigma_vm');
end

plot_summary_ratio(output_root);
plot_summary_solver(output_root);
plot_summary_time(output_root);

% 补充图像：残差曲线、理论 Kt 对比、网格质量、稀疏结构、内存、载荷校核等。
if exist('plot_extra_figures', 'file') == 2
    plot_extra_figures(output_root);
end

fprintf('全部图像绘制完成。\n');
