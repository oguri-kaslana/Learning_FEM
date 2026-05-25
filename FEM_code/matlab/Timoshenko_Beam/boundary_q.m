function F_global = boundary_q(F_global, P, node_q, values_q)
% ============================================================
% boundary_q
%
% 功能：
%   对 Timoshenko 梁单元施加横向均布载荷，并用 Gauss 积分形成等效节点力。
%
% 输入：
%   F_global : 全局载荷向量
%   P        : 节点坐标矩阵
%   node_q   : 受均布载荷作用的单元节点，n × 2
%              node_q(k,:) = [node1, node2]
%   values_q : 均布载荷大小，n × 1
%
% 输出：
%   F_global : 施加载荷后的全局载荷向量
% ============================================================

    [Gp, Gw] = calculate_Gauss(2);

    number_q = size(node_q, 1);
    number_gauss = length(Gp);

    for k = 1:number_q

        node1 = node_q(k, 1);
        node2 = node_q(k, 2);
        q = values_q(k);

        p1 = P(node1, :);
        p2 = P(node2, :);

        J = create_Jacobi(p1, p2);

        Fe = zeros(4, 1);

        for g = 1:number_gauss

            xi = Gp(g);

            N1 = 0.5 * (1 - xi);
            N2 = 0.5 * (1 + xi);

            % =================================================
            % 与 Euler 梁不同：
            % Euler Hermite 梁的均布载荷等效节点力会产生节点力和等效节点弯矩；
            % 二节点 Timoshenko 梁中，横向位移 w 采用线性插值，
            % 均布横向载荷只直接作用在 w 自由度上。
            % 因此这里使用 Nq = [N1; 0; N2; 0]。
            % =================================================
            Nq = [
                N1;
                0;
                N2;
                0
            ];

            Fe = Fe + Nq * q * J * Gw(g);

        end

        dof = [
            2 * node1 - 1;
            2 * node1;
            2 * node2 - 1;
            2 * node2
        ];

        F_global(dof) = F_global(dof) + Fe;

    end

end
