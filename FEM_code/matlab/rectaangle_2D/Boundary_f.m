function f = Boundary_f( ...
    points,topolocy, ...
    nodes, values)
    number_nodes = size(points,1);
    f = zeros(number_nodes,2);
    for i = 1:length(nodes)
        node_id = nodes(i);
        f(node_id,1) = f(node_id,1) + values(i,1);
        f(node_id,2) = f(node_id,2) + values(i,2);
    end

end
