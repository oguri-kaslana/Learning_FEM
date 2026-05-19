function [k_mod,f_mod] = Boundary_Disp( ...
    type_dimesion,...
    nodes,values, ...
    matrix_k_global,boundary_force)
    if type_dimesion==1
        matrix_k_global(nodes,:) = 0;
        matrix_k_global(:,nodes) = 0;
        matrix_k_global(nodes,nodes) = 1;
        k_mod = matrix_k_global;

        boundary_force(nodes) = values;
        f_mod = boundary_force;
    end
end