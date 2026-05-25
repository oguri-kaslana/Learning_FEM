function F_global = boundary_f(F_global, node_f, values_f)
    for k = 1:size(node_f, 1)

        node_id = node_f(k, 1);
        direction = node_f(k, 2);

        dof = 2 * node_id - 2 + direction;

        F_global(dof) = F_global(dof) + values_f(k);

    end

end