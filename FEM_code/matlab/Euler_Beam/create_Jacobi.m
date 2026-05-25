function J = create_Jacobi(p1,p2)
    x = P2(1,1)-p1(1,1);
    y = P2(1,2)-P1(1,2);
    J = sqer(x^2+y^2);
end