function qb = Boundary_Q( ...
    type_dimesion,element, ...
    points,topolocy,q_fun,q_fun_order)
    if type_dimesion==1
        qb = zeros(size(points,1),1);
        nodes = topolocy(element,:);
        x1 = points(nodes(1));
        x2 = points(nodes(2));
        X = x2 - x1;
        J = X / 2;

        if q_fun_order ==0
           qb(nodes) =  q_fun*X/2*[1;1]; 
        elseif q_fun_order ==1
            [Gauss_weights,Gauss_points] = create_Gauss(1,2);
            qb_value1 = 0;
            qb_value2 = 0;
            for i = 1:length(Gauss_points)
                xi = Gauss_points(i);
                weight = Gauss_weights(i);
                N1 = (1-xi)/2;
                N2 = (1+xi)/2;
                % 等参单元的优势
                x = N1*x1+N2*x2;
                q_value = q_fun(x);
                qb_value1 = qb_value1 + weight*N1*q_value*J;
                qb_value2 = qb_value2 + weight*N2*q_value*J;
            end
            qb(nodes) = [qb_value1;qb_value2];
        end
        
        
    end
end