function [k_mod, f_mod] = Boundary_u( ...
    points, topolocy, ...
    k_global, f, nodes, values)

    number_nodes = size(points,1);
    number_dofs = 2 * number_nodes;

    k_mod = k_global;
    F = zeros(number_dofs,1);

    for i = 1:number_nodes
        F(2*i-1) = f(i,1);
        F(2*i)   = f(i,2);
    end

    for i = 1:length(nodes)

        node_id = nodes(i);

        ux_dof = 2*node_id - 1;
        uy_dof = 2*node_id;

        ux_value = values(i,1);
        uy_value = values(i,2);

        % 处理 x 方向位移约束
        F = F - k_mod(:,ux_dof) * ux_value;

        k_mod(ux_dof,:) = 0;
        k_mod(:,ux_dof) = 0;
        k_mod(ux_dof,ux_dof) = 1;

        F(ux_dof) = ux_value;

        % 处理 y 方向位移约束
        F = F - k_mod(:,uy_dof) * uy_value;

        k_mod(uy_dof,:) = 0;
        k_mod(:,uy_dof) = 0;
        k_mod(uy_dof,uy_dof) = 1;

        F(uy_dof) = uy_value;

    end

    f_mod = zeros(number_nodes,2);

    for i = 1:number_nodes
        f_mod(i,1) = F(2*i-1);
        f_mod(i,2) = F(2*i);
    end

end