function B_local = create_B_local(xi, L)
% ============================================================
% create_B_local
%
% 功能：
%   构造二节点 Timoshenko 梁单元在给定自然坐标 xi 处的弯曲 B 和剪切 B。
%
% 输入：
%   xi : 自然坐标，范围 [-1, 1]
%   L  : 当前单元长度
%
% 输出：
%   B_local : 2 × 4 矩阵
%             第 1 行 Bb：曲率-位移矩阵，kappa = dtheta/dx
%             第 2 行 Bs：剪切应变-位移矩阵，gamma = dw/dx - theta
%
% 自由度顺序：
%   Ue = [w1; theta1; w2; theta2]
% ============================================================

    if L <= 0
        error('Element length L must be positive.');
    end

    % 线性形函数
    N1 = 0.5 * (1 - xi);
    N2 = 0.5 * (1 + xi);

    % 线性形函数对实际坐标 x 的导数
    dN1_dx = -1 / L;
    dN2_dx =  1 / L;

    % ========================================================
    % 与 Euler 梁不同：
    % Euler-Bernoulli 梁中 theta = dw/dx，曲率由 w 的二阶导数给出；
    % Timoshenko 梁中 theta 是独立自由度，曲率为 kappa = dtheta/dx。
    % 因此弯曲 B 只作用在 theta 自由度上。
    % ========================================================
    Bb = [0, dN1_dx, 0, dN2_dx];

    % ========================================================
    % 与 Euler 梁不同：
    % Timoshenko 梁考虑剪切变形，剪切应变定义为
    % gamma = dw/dx - theta。
    % 因此剪切 B 同时包含 w 的一阶导数项和 theta 的形函数项。
    % ========================================================
    Bs = [dN1_dx, -N1, dN2_dx, -N2];

    B_local = [Bb; Bs];

end
