function K_local = create_K_local(E, I, P)

    [Gp, Gw] = calculate_Gauss(2);

    C = create_material(E, I);

    B_global = create_B_global(P);

    number_elements = size(P, 1) - 1;
    number_gauss = length(Gp);

    K_local = zeros(number_elements, 4, 4);

    for e = 1:number_elements

        p1 = P(e, :);
        p2 = P(e + 1, :);

        J = create_Jacobi(p1, p2);

        Ke = zeros(4, 4);

        for g = 1:number_gauss

            B = squeeze(B_global(e, g, :))';

            Ke = Ke + B' * C * B * J * Gw(g);

        end

        K_local(e, :, :) = Ke;

    end

end