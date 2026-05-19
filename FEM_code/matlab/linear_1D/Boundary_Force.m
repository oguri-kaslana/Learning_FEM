function f = Boundary_Force( ...
    type_dimesion,number_elements,...
    nodes,values)
    if type_dimesion==1
        f = zeros(number_elements+1,1);
        f(nodes,1) = values;
    end
end