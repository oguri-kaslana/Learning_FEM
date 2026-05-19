function B= create_B_triangle( ...
    points,topolocy, ...
    jacobi) 
    
    number_elements = size(topolocy,1);
    dN_local = [1,0,-1;0,1,-1];
    B = zeros(3,6,number_elements);
    for i = 1:number_elements
        jacobi_local = jacobi(:,:,i);
        dN_global = jacobi_local\dN_local;
        B(:,:,i) = [
            dN_global(1,1),0,dN_global(1,2),0,dN_global(1,3),0;
            0,dN_global(2,1),0,dN_global(2,2),0,dN_global(2,3);
            dN_global(2,1),dN_global(1,1),dN_global(2,2),dN_global(1,2),dN_global(2,3),dN_global(1,3)]; 
    end

end