function kappa = calculate_kappa(P, T, U_global)
% ============================================================
% calculate_kappa
%
% 功能：
%   计算 Timoshenko 梁每个单元在弯曲 Gauss 点处的曲率 kappa。
% ============================================================

    [Bb_global, ~] = create_B_global(P);

    number_elements = size(T, 1);
    number_gauss_Bb = size(Bb_global, 2);

    kappa = zeros(number_elements, number_gauss_Bb);

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

        for g = 1:number_gauss_Bb

            Bb = reshape(Bb_global(e, g, :), 1, 4);

            % =================================================
            % 与 Euler 梁不同：
            % Euler 梁曲率为 kappa = d²w/dx²；
            % Timoshenko 梁曲率为 kappa = dtheta/dx。
            % =================================================
            kappa(e, g) = Bb * Ue;

        end

    end

end
