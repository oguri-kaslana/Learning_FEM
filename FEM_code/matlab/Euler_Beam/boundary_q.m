function F_global = boundary_q(F_global, P, node_q, values_q)

    [Gp, Gw] = calculate_Gauss(2);

    number_q = size(node_q, 1);
    number_gauss = length(Gp);

    for k = 1:number_q

        node1 = node_q(k, 1);
        node2 = node_q(k, 2);

        q = values_q(k);

        p1 = P(node1, :);
        p2 = P(node2, :);

        J = create_Jacobi(p1, p2);

        L = 2 * J;

        Fe = zeros(4, 1);

        for g = 1:number_gauss

            xi = Gp(g);

            N = [
                1/4 * (2 - 3*xi + xi^3);
                L/8 * (1 - xi - xi^2 + xi^3);
                1/4 * (2 + 3*xi - xi^3);
                L/8 * (-1 - xi + xi^2 + xi^3)
            ];

            Fe = Fe + N * q * J * Gw(g);

        end

        dof = [
            2 * node1 - 1;
            2 * node1;
            2 * node2 - 1;
            2 * node2
        ];

        F_global(dof) = F_global(dof) + Fe;

    end

end