%% main.m
% ============================================================
% Timoshenko 梁单元有限元程序调试主文件
%
% 模型：左端固定悬臂梁，右端或全梁施加载荷
% 自由度：每个节点 [w, theta]
%       w     : 横向挠度
%       theta : 独立截面转角
%
% 与 Euler-Bernoulli 梁不同：
%   1. Timoshenko 梁考虑剪切变形；
%   2. theta 不再强制等于 dw/dx；
%   3. 曲率 kappa = dtheta/dx；
%   4. 剪切应变 gamma = dw/dx - theta；
%   5. 刚度包含弯曲刚度 EI 和剪切刚度 kappa_s*G*A。
%
% 本主程序测试三种工况：
%   1. point_force  : 右端集中力
%   2. end_moment   : 右端集中弯矩
%   3. uniform_load : 全梁均布载荷
% ============================================================

clear;
clc;
close all;

%% 1. 基本参数设置
L_total = 1.0;          % 梁总长度，单位 m
b = 0.02;               % 矩形截面宽度，单位 m
h = 0.08;               % 矩形截面高度，单位 m

E = 2.10e11;            % 弹性模量，单位 Pa
nu = 0.30;              % 泊松比
G = E / (2 * (1 + nu)); % 剪切模量，单位 Pa

A = b * h;              % 截面面积，单位 m^2
I = b * h^3 / 12;       % 矩形截面对中性轴的惯性矩，单位 m^4
kappa_s = 5 / 6;        % 矩形截面剪切修正系数

seeds = 0.05;           % 单元长度，单位 m
slender_ratio = L_total / h;

% 载荷大小。负号表示沿负 w 方向。
P0 = -1000.0;           % 右端集中力，单位 N
M0 = -1000.0;           % 右端集中弯矩，单位 N*m
q0 = -1000.0;           % 全梁均布载荷，单位 N/m

fprintf('Timoshenko beam debug model\n');
fprintf('L/h = %.4f\n', slender_ratio);
fprintf('E = %.6e Pa, G = %.6e Pa\n', E, G);
fprintf('A = %.6e m^2, I = %.6e m^4, kappa_s = %.6f\n', A, I, kappa_s);

%% 2. 创建节点坐标矩阵 P 和单元拓扑矩阵 T
[P, T] = create_PT(L_total, seeds);

number_nodes    = size(P, 1);
number_elements = size(T, 1);
total_dofs      = 2 * number_nodes;

left_node  = 1;
right_node = number_nodes;

%% 3. 定义测试工况
case_list = {
    'point_force';
    'end_moment';
    'uniform_load'
};

results = struct();

%% 4. 循环计算三种工况
for case_id = 1:length(case_list)

    case_name = case_list{case_id};

    fprintf('\n============================================================\n');
    fprintf('当前计算工况：%s\n', case_name);
    fprintf('============================================================\n');

    %% 4.1 创建整体刚度矩阵和整体载荷向量
    K_global = create_K_global(P, T, E, I, kappa_s, G, A);
    F_global = zeros(total_dofs, 1);

    %% 4.2 施加载荷边界条件
    switch case_name

        case 'point_force'
            node_f = [right_node, 1];
            values_f = P0;
            F_global = boundary_f(F_global, node_f, values_f);

        case 'end_moment'
            node_f = [right_node, 2];
            values_f = M0;
            F_global = boundary_f(F_global, node_f, values_f);

        case 'uniform_load'
            node_q = T;
            values_q = q0 * ones(number_elements, 1);
            F_global = boundary_q(F_global, P, node_q, values_q);

        otherwise
            error('Unknown load case.');

    end

    %% 4.3 保存原始 K 和 F，用于计算支反力
    K_global_original = K_global;
    F_global_original = F_global;

    %% 4.4 左端固定：w = 0, theta = 0
    node_u = [
        left_node, 1;
        left_node, 2
    ];

    values_u = [0; 0];

    [K_global_bc, F_global_bc] = boundary_u(K_global, F_global, node_u, values_u);

    %% 4.5 求解全局位移
    U_global = calculate_u(K_global_bc, F_global_bc);

    %% 4.6 后处理：支反力、转角、曲率、弯矩、剪切应变、剪力
    RF_global = calculate_RF(K_global_original, F_global_original, U_global);
    angle     = calculate_angle(U_global);
    kappa     = calculate_kappa(P, T, U_global);
    M         = calculate_M(E, I, kappa);

    % 与 Euler 梁不同：Timoshenko 梁需要额外计算剪切应变 gamma 和剪力 Q。
    gamma     = calculate_gamma(P, T, U_global);
    Q         = calculate_Q(kappa_s, G, A, gamma);

    w = U_global(1:2:end);

    %% 4.7 Timoshenko 梁解析解对比：右端挠度与右端转角
    switch case_name

        case 'point_force'
            % Timoshenko 端部集中力悬臂梁：弯曲挠度 + 剪切挠度
            w_right_exact     = P0 * L_total^3 / (3 * E * I) + P0 * L_total / (kappa_s * G * A);
            theta_right_exact = P0 * L_total^2 / (2 * E * I);

        case 'end_moment'
            % 端部纯弯矩下剪力为 0，因此无剪切挠度项
            w_right_exact     = M0 * L_total^2 / (2 * E * I);
            theta_right_exact = M0 * L_total / (E * I);

        case 'uniform_load'
            % Timoshenko 均布载荷悬臂梁：弯曲挠度 + 剪切挠度
            w_right_exact     = q0 * L_total^4 / (8 * E * I) + q0 * L_total^2 / (2 * kappa_s * G * A);
            theta_right_exact = q0 * L_total^3 / (6 * E * I);

    end

    w_right_num     = w(right_node);
    theta_right_num = angle(right_node);

    %% 4.8 输出调试信息
    fprintf('节点数 number_nodes       = %d\n', number_nodes);
    fprintf('单元数 number_elements    = %d\n', number_elements);
    fprintf('总自由度 total_dofs       = %d\n', total_dofs);

    fprintf('\n右端挠度 w_L：\n');
    fprintf('  数值解 = %.12e\n', w_right_num);
    fprintf('  解析解 = %.12e\n', w_right_exact);
    fprintf('  误差   = %.12e\n', w_right_num - w_right_exact);

    fprintf('\n右端转角 theta_L：\n');
    fprintf('  数值解 = %.12e\n', theta_right_num);
    fprintf('  解析解 = %.12e\n', theta_right_exact);
    fprintf('  误差   = %.12e\n', theta_right_num - theta_right_exact);

    fprintf('\n左端支反力与反力矩：\n');
    fprintf('  RF_w_left     = %.12e\n', RF_global(1));
    fprintf('  RF_theta_left = %.12e\n', RF_global(2));

    fprintf('\n各单元弯曲 Gauss 点曲率 kappa：\n');
    disp(kappa);

    fprintf('各单元弯曲 Gauss 点弯矩 M：\n');
    disp(M);

    fprintf('各单元剪切 Gauss 点剪切应变 gamma：\n');
    disp(gamma);

    fprintf('各单元剪切 Gauss 点剪力 Q：\n');
    disp(Q);


end

