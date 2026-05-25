function gamma = calculate_gamma(P, T, U_global)
% ============================================================
% calculate_gamma
%
% 功能：
%   计算 Timoshenko 梁每个单元在剪切 Gauss 点处的剪切应变 gamma。
%
% 说明：
%   Euler-Bernoulli 梁忽略剪切变形，不需要 gamma；
%   Timoshenko 梁需要计算 gamma = dw/dx - theta。
% ============================================================

    [~, Bs_global] = create_B_global(P);

    number_elements = size(T, 1);
    number_gauss_Bs = size(Bs_global, 2);

    gamma = zeros(number_elements, number_gauss_Bs);

    for e = 1:number_elements

        node1 = T(e, 1);
        node2 = T(e, 2);

        dof = [
            2 * node1 - 1;
            2 * node1;
            2 * node2 - 1;
            2 * node2
        ];

        Ue = U_global(dof);

        for g = 1:number_gauss_Bs
            Bs = reshape(Bs_global(e, g, :), 1, 4);
            gamma(e, g) = Bs * Ue;
        end

    end

end
