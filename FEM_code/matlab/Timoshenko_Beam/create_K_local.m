function K_local = create_K_local(E, I, kappa_s, G, A, P)
% ============================================================
% create_K_local
%
% 功能：
%   计算所有 Timoshenko 梁单元的局部刚度矩阵。
%
% 输入：
%   E       : 弹性模量
%   I       : 截面惯性矩
%   kappa_s : 剪切修正系数，矩形截面常取 5/6
%   G       : 剪切模量
%   A       : 截面面积
%   P       : 节点坐标矩阵
%
% 输出：
%   K_local : number_elements × 4 × 4
% ============================================================

    [~, Gw_Bb] = calculate_Gauss(2);
    [~, Gw_Bs] = calculate_Gauss(1);

    C = create_material(E, I, kappa_s, G, A);

    [Bb_global, Bs_global] = create_B_global(P);

    number_elements = size(P, 1) - 1;
    number_gauss_Bb = size(Bb_global, 2);
    number_gauss_Bs = size(Bs_global, 2);

    K_local = zeros(number_elements, 4, 4);

    for e = 1:number_elements

        p1 = P(e, :);
        p2 = P(e + 1, :);
        J = create_Jacobi(p1, p2);

        Ke = zeros(4, 4);

        % ====================================================
        % 与 Euler 梁不同：
        % Timoshenko 梁刚度分为弯曲刚度和剪切刚度两部分。
        % 这里先计算弯曲部分：Ke_b = ∫ Bb' * EI * Bb dx。
        % ====================================================
        for g = 1:number_gauss_Bb
            Bb = reshape(Bb_global(e, g, :), 1, 4);
            Ke = Ke + Bb' * C(1, 1) * Bb * J * Gw_Bb(g);
        end

        % ====================================================
        % 与 Euler 梁不同：
        % Timoshenko 梁还需要剪切刚度：Ke_s = ∫ Bs' * kGA * Bs dx。
        % 这里采用 1 点 Gauss 减缩积分，以减轻剪切锁死。
        % ====================================================
        for g = 1:number_gauss_Bs
            Bs = reshape(Bs_global(e, g, :), 1, 4);
            Ke = Ke + Bs' * C(2, 2) * Bs * J * Gw_Bs(g);
        end

        K_local(e, :, :) = Ke;

    end

end
