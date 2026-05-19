function k_global = create_K_global(topolocy, points, k_local)

    number_dofs = 2 * size(points, 1);
    k_global = zeros(number_dofs, number_dofs);

    for i = 1:size(topolocy, 1)

        nodes = topolocy(i, :);

        dofs = [
            2*nodes(1)-1, 2*nodes(1), ...
            2*nodes(2)-1, 2*nodes(2), ...
            2*nodes(3)-1, 2*nodes(3)
        ];

        k_global(dofs, dofs) = k_global(dofs, dofs) + k_local(:,:,i);

    end

end