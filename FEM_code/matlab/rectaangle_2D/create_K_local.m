function k_local = create_K_local(topolocy,points, ...
    B,c,det_jacobi)
    number_elements = size(topolocy,1);

    k_local = zeros(6,6,number_elements);


    for i = 1:number_elements

        Ke = 1/2*B(:,:,i)'*c*B(:,:,i)*det_jacobi(i);
        
        k_local(:,:,i) = Ke;

    end

end