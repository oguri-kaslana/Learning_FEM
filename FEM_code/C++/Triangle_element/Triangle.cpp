#pragma once
#include <vector>
#include <map>
#include <utility>
#include <string>
#include "Calcaulate.h"
#include <cmath>
using namespace std;

//生成模型的节点的整体编号
vector<vector<int>> create_nodes(double L, double H, double seeds)
{
    int nx = int((L) / (seeds)) + 1;
    int ny = int((H) / seeds) + 1;
    vector<vector<int>> nodes(ny, vector<int>(nx));
    for (int i = 0; i < ny; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            nodes[ny - i - 1][j] = j + 1 + nx * i;
        }
    }
    return nodes;
}

//生成节点坐标阵P
vector<vector<double>> create_Points(double L, double H, double seeds)
{
    int nx = int((L) / (seeds)) + 1;
    double dx = L / (nx - 1);
    int ny = int((H) / seeds) + 1;
    double dy = H / (ny - 1);
    int number_nodes = nx * ny;

    vector<vector<double>> P(number_nodes, vector<double>(2));
    for (int i = 0; i < ny; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            P[i * nx + j][0] = j * dx;
            P[i * nx + j][1] = i * dy;
        }
    }

    return P;
}

//生成连接矩阵T
vector<vector<int>> create_Topolocy(vector<vector<int>> nodes)
{
    int ny = nodes.size();
    int nx = nodes[0].size();

    int number_elements = 2 * (nx - 1) * (ny - 1);

    vector<vector<int>> T(number_elements, vector<int>(3));

    int e = 0;

    // 从模型左下角开始，逐行向右，再向上生成单元
    for (int i = ny - 1; i > 0; i--)
    {
        for (int j = 0; j < nx - 1; j++)
        {
            int n1 = nodes[i][j];         // 左下
            int n2 = nodes[i][j + 1];     // 右下
            int n3 = nodes[i - 1][j + 1]; // 右上
            int n4 = nodes[i - 1][j];     // 左上

            // 第一个三角形
            T[e][0] = n1;
            T[e][1] = n2;
            T[e][2] = n3;
            e++;

            // 第二个三角形
            T[e][0] = n1;
            T[e][1] = n3;
            T[e][2] = n4;
            e++;
        }
    }

    return T;
}

//生成Jacobi
vector<vector<vector<double>>> create_Jacobi(vector<vector<double>> P,
    vector<vector<int>> T)
{
    int number_elements = T.size();

    vector<vector<vector<double>>> Jacobi(
        2, vector<vector<double>>(2, vector<double>(number_elements))
    );

    for (int e = 0; e < number_elements; e++)
    {
        int n1 = T[e][0] - 1;
        int n2 = T[e][1] - 1;
        int n3 = T[e][2] - 1;

        double x1 = P[n1][0];
        double y1 = P[n1][1];

        double x2 = P[n2][0];
        double y2 = P[n2][1];

        double x3 = P[n3][0];
        double y3 = P[n3][1];

        Jacobi[0][0][e] = x2 - x1;
        Jacobi[0][1][e] = x3 - x1;
        Jacobi[1][0][e] = y2 - y1;
        Jacobi[1][1][e] = y3 - y1;
    }

    return Jacobi;
}

//生成等参单元的形函数的导数
vector<vector<double>> create_dN_local()
{
    // 三节点线性三角形等参单元
    // N1 = 1 - ksi - eta
    // N2 = ksi
    // N3 = eta

    vector<vector<double>> dN_local(2, vector<double>(3));

    dN_local[0][0] = -1.0;
    dN_local[0][1] = 1.0;
    dN_local[0][2] = 0.0;

    dN_local[1][0] = -1.0;
    dN_local[1][1] = 0.0;
    dN_local[1][2] = 1.0;

    return dN_local;
}

//组装B矩阵
vector<vector<vector<double>>> create_B_triangle(vector<vector<double>> dN_local,
    vector<vector<vector<double>>> Jacobi)
{
    int number_elements = Jacobi[0][0].size();

    // B[行][列][单元编号]
    vector<vector<vector<double>>> B(
        3, vector<vector<double>>(6, vector<double>(number_elements, 0.0))
    );

    for (int e = 0; e < number_elements; e++)
    {
        // 取出第 e 个单元的 Jacobi 矩阵
        vector<vector<double>> J(2, vector<double>(2, 0.0));

        J[0][0] = Jacobi[0][0][e];
        J[0][1] = Jacobi[0][1][e];
        J[1][0] = Jacobi[1][0][e];
        J[1][1] = Jacobi[1][1][e];

        // invJT = J^(-T)
        vector<vector<double>> invJ = inverse_2x2(J);
        vector<vector<double>> invJT = matrix_transpose(invJ);

        // dN_global = J^(-T) * dN_local
        vector<vector<double>> dN_global = matrix_multiply(invJT, dN_local);

        // 组装B矩阵
        B[0][0][e] = dN_global[0][0];
        B[0][2][e] = dN_global[0][1];
        B[0][4][e] = dN_global[0][2];

        B[1][1][e] = dN_global[1][0];
        B[1][3][e] = dN_global[1][1];
        B[1][5][e] = dN_global[1][2];

        B[2][0][e] = dN_global[1][0];
        B[2][1][e] = dN_global[0][0];

        B[2][2][e] = dN_global[1][1];
        B[2][3][e] = dN_global[0][1];

        B[2][4][e] = dN_global[1][2];
        B[2][5][e] = dN_global[0][2];
    }

    return B;
}

//组装平面应力D矩阵
vector<vector<double>> create_D_plane_stress(double E, double nu)
{
    vector<vector<double>> D(3, vector<double>(3, 0.0));

    double c = E / (1.0 - nu * nu);

    D[0][0] = c;
    D[0][1] = c * nu;

    D[1][0] = c * nu;
    D[1][1] = c;

    D[2][2] = c * (1.0 - nu) / 2.0;

    return D;
}

//组装K_local
vector<vector<vector<double>>> create_K_local(
    vector<vector<vector<double>>> B,
    vector<vector<vector<double>>> Jacobi,
    vector<vector<double>> D,
    double thickness)
{
    int number_elements = B[0][0].size();
    vector<vector<vector<double>>> K_local(
        6, vector<vector<double>>(6, vector<double>(number_elements, 0.0))
    );

    for (int e = 0; e < number_elements; e++)
    {
        // 取出第e个单元的B矩阵
        vector<vector<double>> Be(3, vector<double>(6, 0.0));

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                Be[i][j] = B[i][j][e];
            }
        }

        // 取出第e个单元的Jacobi矩阵
        vector<vector<double>> J(2, vector<double>(2, 0.0));

        J[0][0] = Jacobi[0][0][e];
        J[0][1] = Jacobi[0][1][e];
        J[1][0] = Jacobi[1][0][e];
        J[1][1] = Jacobi[1][1][e];

        double detJ = determinant_2x2(J);


        // K_local = B^T * D * B * thickness * area
        vector<vector<double>> BT = matrix_transpose(Be);
        vector<vector<double>> DB = matrix_multiply(D, Be);
        vector<vector<double>> BTDB = matrix_multiply(BT, DB);

        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                K_local[i][j][e] = BTDB[i][j] * thickness * 0.5 * abs(detJ);
            }
        }
    }

    return K_local;
}

//组装K_global
vector<vector<double>> create_K_global(
    vector<vector<vector<double>>> K,
    vector<vector<int>> T,
    vector<vector<double>> P)
{
    int number_nodes = P.size();
    int number_dof = number_nodes * 2;
    int number_elements = T.size();

    vector<vector<double>> K_global(
        number_dof, vector<double>(number_dof, 0.0)
    );

    for (int e = 0; e < number_elements; e++)
    {
        int n1 = T[e][0] - 1;
        int n2 = T[e][1] - 1;
        int n3 = T[e][2] - 1;

        vector<int> element_dof(6);

        element_dof[0] = 2 * n1;
        element_dof[1] = 2 * n1 + 1;

        element_dof[2] = 2 * n2;
        element_dof[3] = 2 * n2 + 1;

        element_dof[4] = 2 * n3;
        element_dof[5] = 2 * n3 + 1;

        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                K_global[element_dof[i]][element_dof[j]] += K[i][j][e];
            }
        }
    }

    return K_global;
}

//六节点二次单元
pair<vector<vector<double>>, vector<vector<int>>> create_T6_mesh(
    vector<vector<double>> P3,
    vector<vector<int>> T3)
{
    vector<vector<double>> P6 = P3;
    vector<vector<int>> T6(T3.size(), vector<int>(6, 0));

    map<pair<int, int>, int> edge_mid_node;

    for (int e = 0; e < T3.size(); e++)
    {
        int n1 = T3[e][0];
        int n2 = T3[e][1];
        int n3 = T3[e][2];

        vector<pair<int, int>> edges = {
            {n1, n2},
            {n2, n3},
            {n3, n1}
        };

        vector<int> mid_nodes(3, 0);

        for (int k = 0; k < 3; k++)
        {
            int a = edges[k].first;
            int b = edges[k].second;

            pair<int, int> edge_key = { min(a, b), max(a, b) };

            if (edge_mid_node.find(edge_key) == edge_mid_node.end())
            {
                double x_mid = 0.5 * (P3[a - 1][0] + P3[b - 1][0]);
                double y_mid = 0.5 * (P3[a - 1][1] + P3[b - 1][1]);

                P6.push_back({ x_mid, y_mid });

                int mid_id = P6.size();
                edge_mid_node[edge_key] = mid_id;
            }

            mid_nodes[k] = edge_mid_node[edge_key];
        }

        T6[e][0] = n1;
        T6[e][1] = n2;
        T6[e][2] = n3;
        T6[e][3] = mid_nodes[0];
        T6[e][4] = mid_nodes[1];
        T6[e][5] = mid_nodes[2];
    }

    return { P6, T6 };
}

vector<vector<double>> create_dN_local_T6(double ksi, double eta)
{
    vector<vector<double>> dN_local(2, vector<double>(6, 0.0));

    double L1 = 1.0 - ksi - eta;
    double L2 = ksi;
    double L3 = eta;

    dN_local[0][0] = 1.0 - 4.0 * L1;
    dN_local[0][1] = 4.0 * L2 - 1.0;
    dN_local[0][2] = 0.0;
    dN_local[0][3] = 4.0 * (L1 - L2);
    dN_local[0][4] = 4.0 * L3;
    dN_local[0][5] = -4.0 * L3;

    dN_local[1][0] = 1.0 - 4.0 * L1;
    dN_local[1][1] = 0.0;
    dN_local[1][2] = 4.0 * L3 - 1.0;
    dN_local[1][3] = -4.0 * L2;
    dN_local[1][4] = 4.0 * L2;
    dN_local[1][5] = 4.0 * (L1 - L3);

    return dN_local;
}

vector<vector<double>> create_B_triangle_T6_single(
    vector<vector<double>> P,
    vector<int> element_nodes,
    double ksi,
    double eta)
{
    int n1 = element_nodes[0] - 1;
    int n2 = element_nodes[1] - 1;
    int n3 = element_nodes[2] - 1;

    vector<vector<double>> J(2, vector<double>(2, 0.0));

    J[0][0] = P[n2][0] - P[n1][0];
    J[0][1] = P[n3][0] - P[n1][0];
    J[1][0] = P[n2][1] - P[n1][1];
    J[1][1] = P[n3][1] - P[n1][1];

    vector<vector<double>> invJ = inverse_2x2(J);
    vector<vector<double>> invJT = matrix_transpose(invJ);

    vector<vector<double>> dN_local = create_dN_local_T6(ksi, eta);

    vector<vector<double>> dN_global = matrix_multiply(invJT, dN_local);

    vector<vector<double>> B(3, vector<double>(12, 0.0));

    for (int i = 0; i < 6; i++)
    {
        B[0][2 * i] = dN_global[0][i];
        B[1][2 * i + 1] = dN_global[1][i];

        B[2][2 * i] = dN_global[1][i];
        B[2][2 * i + 1] = dN_global[0][i];
    }

    return B;
}

vector<vector<vector<double>>> create_K_local_T6(
    vector<vector<double>> P,
    vector<vector<int>> T6,
    vector<vector<double>> D,
    double thickness)
{
    int number_elements = T6.size();

    vector<vector<vector<double>>> K_local(
        12, vector<vector<double>>(12, vector<double>(number_elements, 0.0))
    );

    vector<vector<double>> gauss_points = {
        {1.0 / 6.0, 1.0 / 6.0},
        {2.0 / 3.0, 1.0 / 6.0},
        {1.0 / 6.0, 2.0 / 3.0}
    };

    vector<double> gauss_weights = {
        1.0 / 6.0,
        1.0 / 6.0,
        1.0 / 6.0
    };

    for (int e = 0; e < number_elements; e++)
    {
        int n1 = T6[e][0] - 1;
        int n2 = T6[e][1] - 1;
        int n3 = T6[e][2] - 1;

        vector<vector<double>> J(2, vector<double>(2, 0.0));

        J[0][0] = P[n2][0] - P[n1][0];
        J[0][1] = P[n3][0] - P[n1][0];
        J[1][0] = P[n2][1] - P[n1][1];
        J[1][1] = P[n3][1] - P[n1][1];

        double detJ = abs(determinant_2x2(J));

        for (int gp = 0; gp < gauss_points.size(); gp++)
        {
            double ksi = gauss_points[gp][0];
            double eta = gauss_points[gp][1];
            double weight = gauss_weights[gp];

            vector<vector<double>> B =
                create_B_triangle_T6_single(P, T6[e], ksi, eta);

            vector<vector<double>> BT = matrix_transpose(B);
            vector<vector<double>> DB = matrix_multiply(D, B);
            vector<vector<double>> BTDB = matrix_multiply(BT, DB);

            for (int i = 0; i < 12; i++)
            {
                for (int j = 0; j < 12; j++)
                {
                    K_local[i][j][e] +=
                        BTDB[i][j] * thickness * detJ * weight;
                }
            }
        }
    }

    return K_local;
}

vector<vector<double>> create_K_global_general(
    vector<vector<vector<double>>> K_local,
    vector<vector<int>> T,
    vector<vector<double>> P)
{
    int number_nodes = P.size();
    int number_dof = number_nodes * 2;
    int number_elements = T.size();
    int number_element_nodes = T[0].size();
    int number_element_dof = number_element_nodes * 2;

    vector<vector<double>> K_global(
        number_dof, vector<double>(number_dof, 0.0)
    );

    for (int e = 0; e < number_elements; e++)
    {
        vector<int> element_dof(number_element_dof, 0);

        for (int i = 0; i < number_element_nodes; i++)
        {
            int node_id = T[e][i] - 1;

            element_dof[2 * i] = 2 * node_id;
            element_dof[2 * i + 1] = 2 * node_id + 1;
        }

        for (int i = 0; i < number_element_dof; i++)
        {
            for (int j = 0; j < number_element_dof; j++)
            {
                K_global[element_dof[i]][element_dof[j]] +=
                    K_local[i][j][e];
            }
        }
    }

    return K_global;
}

vector<double> get_element_displacement_general(
    vector<double> U,
    vector<int> element_nodes)
{
    int number_element_nodes = element_nodes.size();

    vector<double> Ue(number_element_nodes * 2, 0.0);

    for (int i = 0; i < number_element_nodes; i++)
    {
        int node_id = element_nodes[i] - 1;

        Ue[2 * i] = U[2 * node_id];
        Ue[2 * i + 1] = U[2 * node_id + 1];
    }

    return Ue;
}

vector<vector<double>> create_strain_T6(
    vector<vector<double>> P,
    vector<vector<int>> T6,
    vector<double> U)
{
    int number_elements = T6.size();

    vector<vector<double>> strain(number_elements, vector<double>(3, 0.0));

    double ksi = 1.0 / 3.0;
    double eta = 1.0 / 3.0;

    for (int e = 0; e < number_elements; e++)
    {
        vector<vector<double>> B =
            create_B_triangle_T6_single(P, T6[e], ksi, eta);

        vector<double> Ue = get_element_displacement_general(U, T6[e]);

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 12; j++)
            {
                strain[e][i] += B[i][j] * Ue[j];
            }
        }
    }

    return strain;
}
