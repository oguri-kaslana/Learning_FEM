function kappa = calculate_kappa(P, T, U_global)
    B_global = create_B_global(P);

    number_elements = size(T, 1);
    number_gauss = size(B_global, 2);

    kappa = zeros(number_elements, number_gauss);

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

        for g = 1:number_gauss

            B = reshape(B_global(e, g, :), 1, 4);

            kappa(e, g) = B * Ue;

        end

    end

end