function K_local = create_K_local(E,I,P)
    [~,Gw] = calculate_Gauss(2);
    C = create_material(E,I);
    B_global = create_B_global(P);
    number_elements = size(P,1)-1;
    K_local = zeros(number_elements,4,4);
    for e = 1:number_elements
        p1 = P(e, :);
        p2 = P(e + 1, :);
        J = create_Jacobi(p1, p2);
        K_local(e,:,:) = C*det(J)*Gw(1)*B_global(e,1,:)*Gw(2)*B_global(e,2,:);
    end


end