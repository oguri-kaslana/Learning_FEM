function [K_global, F_global] = boundary_u(K_global, F_global, node_u, values_u)
    number_dofs = size(K_global, 1);

    for k = 1:size(node_u, 1)

        node_id = node_u(k, 1);
        direction = node_u(k, 2);

        dof = 2 * node_id - 2 + direction;

        value = values_u(k);

        for i = 1:number_dofs
            F_global(i) = F_global(i) - K_global(i, dof) * value;
        end

        K_global(dof, :) = 0;
        K_global(:, dof) = 0;

        K_global(dof, dof) = 1;

        F_global(dof) = value;

    end

end