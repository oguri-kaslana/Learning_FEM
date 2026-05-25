%% main.m
% ============================================================
% Euler-Bernoulli 梁单元有限元程序调试主文件
%
% 模型：左端固定悬臂梁，右端或全梁施加载荷
% 自由度：每个节点 [w, theta]
%       w     : 横向挠度
%       theta : 截面转角
%
% 本主程序分别测试三种工况：
%   1. point_force  : 右端集中力
%   2. end_moment   : 右端集中弯矩
%   3. uniform_load : 全梁均布载荷，用 boundary_q 通过高斯积分形成等效节点力
%
% 符号约定：
%   w 正方向、theta 正方向与单元插值和刚度矩阵保持一致。
%   下面 P0、M0、q0 取负值表示反方向载荷，便于模拟向下荷载。
% ============================================================

clear;
clc;
close all;

%% 1. 基本参数设置
% ------------------------------------------------------------
% 为便于调试，先取 E = 1, I = 1, L = 1。
% 这样解析解数值简单，方便检查程序是否正确。
% 程序跑通后，可以替换成真实工程参数。
% ------------------------------------------------------------
L_total = 1.0;              % 梁总长度，单位 m
b = 0.02;             % 截面宽度，单位 m
h = 0.02;             % 截面高度，单位 m

% Euler-Bernoulli 梁要求 L/h 足够大
slender_ratio = L_total / h;

E = 2.10e11;          % 弹性模量，单位 Pa，钢材典型值
I = b * h^3 / 12;     % 矩形截面对中性轴的惯性矩，单位 m^4

seeds = 0.05;         % 单元长度，单位 m

% 三种载荷大小
P0 = -1.0;              % 右端集中力，作用在 w 自由度
M0 = -1.0;              % 右端集中弯矩，作用在 theta 自由度
q0 = -1.0;              % 全梁均布载荷，单位长度载荷

%% 2. 创建节点坐标矩阵 P 和单元拓扑矩阵 T
[P, T] = create_PT(L_total, seeds);

number_nodes    = size(P, 1);
number_elements = size(T, 1);
total_dofs      = 2 * number_nodes;

left_node  = 1;
right_node = number_nodes;

%% 3. 定义需要测试的工况
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
    K_global = create_K_global(P, T, E, I);
    F_global = zeros(total_dofs, 1);

    %% 4.2 施加载荷边界条件
    switch case_name

        case 'point_force'
            % ------------------------------------------------
            % 右端集中力：施加在右端节点的 w 自由度
            % node_f(k,:) = [node_id, direction]
            % direction = 1 表示 w 自由度
            % direction = 2 表示 theta 自由度
            % ------------------------------------------------
            node_f = [right_node, 1];
            values_f = P0;
            F_global = boundary_f(F_global, node_f, values_f);

        case 'end_moment'
            % ------------------------------------------------
            % 右端集中弯矩：施加在右端节点的 theta 自由度
            % ------------------------------------------------
            node_f = [right_node, 2];
            values_f = M0;
            F_global = boundary_f(F_global, node_f, values_f);

        case 'uniform_load'
            % ------------------------------------------------
            % 全梁均布载荷：对每一个梁单元施加 q0
            % boundary_q 内部使用高斯积分计算等效节点力
            % node_q(k,:) = [node1, node2]
            % ------------------------------------------------
            node_q = T;
            values_q = q0 * ones(number_elements, 1);
            F_global = boundary_q(F_global, P, node_q, values_q);

        otherwise
            error('未知工况类型。');

    end

    %% 4.3 保存原始 K 和 F，用于后续计算支反力
    K_global_original = K_global;
    F_global_original = F_global;

    %% 4.4 施加位移边界条件：左端固定
    % 左端固定表示：
    %   w_1     = 0
    %   theta_1 = 0
    node_u = [
        left_node, 1;
        left_node, 2
    ];

    values_u = [
        0;
        0
    ];

    [K_global_bc, F_global_bc] = boundary_u(K_global, F_global, node_u, values_u);

    %% 4.5 求解全局位移向量
    U_global = calculate_u(K_global_bc, F_global_bc);

    %% 4.6 计算支反力、转角、曲率、弯矩
    RF_global = calculate_RF(K_global_original, F_global_original, U_global);
    angle     = calculate_angle(U_global);
    kappa     = calculate_kappa(P, T, U_global);
    M         = calculate_M(E, I, kappa);

    % 提取节点挠度
    w = U_global(1:2:end);

    %% 4.7 解析解对比：只检查右端挠度和右端转角
    switch case_name

        case 'point_force'
            w_right_exact     = P0 * L_total^3 / (3 * E * I);
            theta_right_exact = P0 * L_total^2 / (2 * E * I);

        case 'end_moment'
            w_right_exact     = M0 * L_total^2 / (2 * E * I);
            theta_right_exact = M0 * L_total / (E * I);

        case 'uniform_load'
            w_right_exact     = q0 * L_total^4 / (8 * E * I);
            theta_right_exact = q0 * L_total^3 / (6 * E * I);

    end

    w_right_num     = w(right_node);
    theta_right_num = angle(right_node);

    %% 4.8 输出主要调试信息
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

    fprintf('\n各单元各高斯点曲率 kappa：\n');
    disp(kappa);

    fprintf('各单元各高斯点弯矩 M：\n');
    disp(M);

 

end
