function [points, topolocy,nodes]  = create_mesh_triangle( ...
    left, right, down, top)
    x = left:1:right;
    y = down:1:top;
    row = length(y);
    col = length(x);
    number_nodes = length(x)*length(y);
    number_elements = 2*(length(x)-1)*(length(y)-1);
    row_elements = 2*(length(x)-1);
    points = zeros(number_nodes,2);
    topolocy = zeros(number_elements,3);
    nodes = zeros(length(y),length(x));
    for i = 1:row
        nodes(row-i+1,:) = 1 + (i-1)*col:1:col + (i-1)*col;
        for j = 1:col
            points(j+(i-1)*col,:) = [x(j),y(i)];
        end
    end
    for i = 1:row-1
        for j = 1:col-1
            topolocy(2*j-1+(i-1)*row_elements,:) = [nodes(row-i+1,j),nodes(row-i+1,j+1),nodes(row-i,j+1)];
            topolocy(2*j+(i-1)*row_elements,:) = [nodes(row-i+1,j),nodes(row-i,j+1),nodes(row-i,j)];
        end
    end
end