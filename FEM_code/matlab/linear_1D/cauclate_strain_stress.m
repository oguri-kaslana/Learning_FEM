function [strain,stress] = cauclate_strain_stress( ...
    number_dimesion,points,topolocy, ...
    u,c)
    if number_dimesion ==1
        number_nodes = size(points,1);
        strain = zeros(number_nodes,1);
        for i = 1:number_nodes-1
            X = points(topolocy(i,:),:);
            l = X(2)-X(1);
            delt_u = u(topolocy(i,2),:)-u(topolocy(i,1),:);
            strain(i,:) = strain(i,:) + delt_u/l;

        end
        stress = strain*c;
    end
end