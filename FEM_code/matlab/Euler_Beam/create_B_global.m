function B_global = create_B_global(P)

    number_elements = size(P, 1) - 1;
    [Gp, Gw] = calculate_Gauss(2);
    number_gauss = length(Gp);
    B_global = zeros(number_elements, number_gauss, 4);
    for e = 1:number_elements
        p1 = P(e, :);
        p2 = P(e + 1, :);
        J = create_Jacobi(p1, p2);
        L = 2*J;
        for g = 1:number_gauss
            xi = Gp(g);
            B_local = create_B_local(xi,L);
            B_global(e, g, :) = B_local / J^2;

        end

    end

end