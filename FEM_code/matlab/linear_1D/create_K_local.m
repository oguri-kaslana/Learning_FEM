function k_local = create_K_local( ...
    number_dimesion,...
    Jacobi,c)
    if number_dimesion==1
        number_elements = size(Jacobi,2);
        for i = 1:number_elements
            J = Jacobi(i,:);
            k_local = [1,-1;-1,1]/(2*J)*c;
        end
    end

end