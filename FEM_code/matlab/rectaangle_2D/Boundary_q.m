function q = Boundary_q( ...
    points,topolocy, thickness,...
    nodes,values) 
    number_nodes = size(points,1);
    q = zeros(number_nodes,2);
    x = points(nodes,1);
    y = points(nodes,2);
    for i = 1: length(nodes)-1
        L = abs(y(i)-y(i+1));
        q(nodes(i),:) =  q(nodes(i),:)+1/2*L*thickness*[values(i,1),values(i,2)];
        q(nodes(i+1),:) = q(nodes(i),:)+1/2*L*thickness*[values(i,1),values(i,2)];
    end
end