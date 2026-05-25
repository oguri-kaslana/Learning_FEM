function [Bb_global, Bs_global] = create_B_global(P)
% ============================================================
% create_B_global
%
% 功能：
%   计算所有单元、所有 Gauss 点处的 Timoshenko 梁弯曲 B 和剪切 B。
%
% 输入：
%   P : 节点坐标矩阵
%
% 输出：
%   Bb_global : number_elements × number_gauss_Bb × 4
%               弯曲 B，计算 kappa = Bb * Ue
%
%   Bs_global : number_elements × number_gauss_Bs × 4
%               剪切 B，计算 gamma = Bs * Ue
% ============================================================

    number_elements = size(P, 1) - 1;

    % ========================================================
    % 与 Euler 梁不同：
    % Timoshenko 梁刚度由弯曲项和剪切项组成。
    % 为减轻剪切锁死，通常采用选择性积分：
    %   弯曲项：2 点 Gauss 积分
    %   剪切项：1 点 Gauss 减缩积分
    % ========================================================
    [Gp_Bb, ~] = calculate_Gauss(2);
    [Gp_Bs, ~] = calculate_Gauss(1);

    number_gauss_Bb = length(Gp_Bb);
    number_gauss_Bs = length(Gp_Bs);

    Bb_global = zeros(number_elements, number_gauss_Bb, 4);
    Bs_global = zeros(number_elements, number_gauss_Bs, 4);

    for e = 1:number_elements

        p1 = P(e, :);
        p2 = P(e + 1, :);

        J = create_Jacobi(p1, p2);
        L = 2 * J;

        for g = 1:number_gauss_Bb
            xi = Gp_Bb(g);
            B_local = create_B_local(xi, L);
            Bb_global(e, g, :) = B_local(1, :);
        end

        for g = 1:number_gauss_Bs
            xi = Gp_Bs(g);
            B_local = create_B_local(xi, L);
            Bs_global(e, g, :) = B_local(2, :);
        end

    end

end
