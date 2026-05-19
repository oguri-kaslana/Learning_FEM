function u = cauclate_U(k_mod, f_mod)
    number_nodes = size(f_mod,1);
    number_dofs = 2 * number_nodes;
    F = zeros(number_dofs,1);
    for i = 1:number_nodes
        F(2*i-1) = f_mod(i,1); 
        F(2*i)   = f_mod(i,2);   
    end
    U = k_mod \ F;
    u = zeros(number_nodes,2);
    for i = 1:number_nodes
        u(i,1) = U(2*i-1);  
        u(i,2) = U(2*i);     
    end

end