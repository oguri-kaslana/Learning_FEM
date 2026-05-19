function k_global = create_K_global( ...
    type_dimesion,points,topolocy,...
    Jacobi,c)
    if type_dimesion==1
        number_elements = size(Jacobi,1);
        number_nodes = size(points,1);
        k_global = zeros(number_nodes,number_nodes);
        for i = 1:number_elements
            matrix_k_local = create_K_local(type_dimesion,Jacobi,c);
            dof = topolocy(i,:);
            k_global(dof,dof) = k_global(dof,dof) + matrix_k_local; 
        end
    end
end