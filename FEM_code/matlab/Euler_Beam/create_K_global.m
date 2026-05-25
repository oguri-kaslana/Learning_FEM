function K_global = create_K_global(P, T, E, I)

    K_local = create_K_local(E, I, P);

    number_nodes = size(P, 1);
    number_elements = size(T, 1);

    total_dofs = 2 * number_nodes;

    K_global = zeros(total_dofs, total_dofs);

    for e = 1:number_elements

        Ke = squeeze(K_local(e, :, :));

        node1 = T(e, 1);
        node2 = T(e, 2);

        dof = [
            2 * node1 - 1, ...
            2 * node1, ...
            2 * node2 - 1, ...
            2 * node2
        ];

        K_global(dof, dof) = K_global(dof, dof) + Ke;

    end

end