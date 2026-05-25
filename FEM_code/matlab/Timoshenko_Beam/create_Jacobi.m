function J = create_Jacobi(p1, p2)

    dx = p2(1) - p1(1);
    dy = p2(2) - p1(2);

    L = sqrt(dx^2 + dy^2);

    J = L / 2;

end