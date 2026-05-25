function [P,T] = create_PT(L,seeds)

    nx = L / seeds;
    ny = 0;

    number_nodes = (nx + 1)*(ny+1);
    number_elements = nx;

    P = zeros(number_nodes, 2);
    T = zeros(number_elements, 2);

    % =============================
    % 生成节点坐标矩阵 P
    % 由节点编号反算网格中的i,j位置
    % =============================
    for node_id = 1:number_nodes

        i = mod(node_id - 1, nx + 1);
        j = floor((node_id - 1) / (nx + 1));

        x = i * seeds;
        y = j * seeds;

        P(node_id, :) = [x, y];

    end

    % =============================
    % 生成单元拓扑矩阵 T
    % =============================
    for element_id = 1:number_elements
        % 计算单元的节点编号坐标(j,i)
        i = mod(element_id - 1, nx) + 1;
        j = floor((element_id - 1) / nx) + 1;

        node1 = (j - 1) * (nx + 1) + i;
        node2 = node1 + 1;

        T(element_id, :) = [node1, node2];

    end

end