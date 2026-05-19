#include "Boundary.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
using namespace std;

void apply_u_boundary(
    vector<vector<double>>& K_global,
    vector<double>& F_global,
    vector<vector<int>> node_u,
    vector<double> values_u)
{
    int number_dof = K_global.size();

    for (int k = 0; k < node_u.size(); k++)
    {
        int node_id = node_u[k][0];
        int direction = node_u[k][1];

        int dof = 2 * (node_id - 1) + (direction - 1);
        double value = values_u[k];

        for (int i = 0; i < number_dof; i++)
        {
            F_global[i] -= K_global[i][dof] * value;
        }

        for (int j = 0; j < number_dof; j++)
        {
            K_global[dof][j] = 0.0;
        }

        for (int i = 0; i < number_dof; i++)
        {
            K_global[i][dof] = 0.0;
        }

        K_global[dof][dof] = 1.0;
        F_global[dof] = value;
    }
}

void apply_f_boundary(
    vector<double>& F_global,
    vector<vector<int>> node_f,
    vector<double> values_f)
{
    for (int k = 0; k < node_f.size(); k++)
    {
        int node_id = node_f[k][0];
        int direction = node_f[k][1];

        int dof = 2 * (node_id - 1) + (direction - 1);

        F_global[dof] += values_f[k];
    }
}

void apply_q_boundary(
    vector<double>& F_global,
    vector<vector<double>> P,
    vector<vector<int>> node_q,
    vector<vector<double>> values_q,
    double thickness)
{
    for (int k = 0; k < node_q.size(); k++)
    {
        int node1 = node_q[k][0];
        int node2 = node_q[k][1];

        double qx = values_q[k][0];
        double qy = values_q[k][1];

        double x1 = P[node1 - 1][0];
        double y1 = P[node1 - 1][1];

        double x2 = P[node2 - 1][0];
        double y2 = P[node2 - 1][1];

        double length = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

        int dof1_x = 2 * (node1 - 1);
        int dof1_y = 2 * (node1 - 1) + 1;

        int dof2_x = 2 * (node2 - 1);
        int dof2_y = 2 * (node2 - 1) + 1;

        F_global[dof1_x] += qx * length * thickness / 2.0;
        F_global[dof1_y] += qy * length * thickness / 2.0;

        F_global[dof2_x] += qx * length * thickness / 2.0;
        F_global[dof2_y] += qy * length * thickness / 2.0;
    }
}

// 自动寻找左端固定节点
vector<vector<int>> find_left_fixed_nodes(vector<vector<double>> P)
{
    vector<vector<int>> node_u;

    double eps = 1.0e-8;

    for (int i = 0; i < P.size(); i++)
    {
        double x = P[i][0];

        if (abs(x - 0.0) < eps)
        {
            int node_id = i + 1;

            node_u.push_back({ node_id, 1 });
            node_u.push_back({ node_id, 2 });
        }
    }

    return node_u;
}

// 自动寻找右端受载边
vector<vector<int>> find_right_boundary_edges(
    vector<vector<double>> P,
    double L)
{
    vector<pair<int, double>> right_nodes;
    vector<vector<int>> node_q;

    double eps = 1.0e-8;

    for (int i = 0; i < P.size(); i++)
    {
        double x = P[i][0];
        double y = P[i][1];

        if (abs(x - L) < eps)
        {
            int node_id = i + 1;
            right_nodes.push_back({ node_id, y });
        }
    }

    sort(right_nodes.begin(), right_nodes.end(),
        [](pair<int, double> a, pair<int, double> b)
        {
            return a.second < b.second;
        });

    if (right_nodes.size() < 2)
    {
        return node_q;
    }

    for (int i = 0; i < right_nodes.size() - 1; i++)
    {
        int node1 = right_nodes[i].first;
        int node2 = right_nodes[i + 1].first;

        node_q.push_back({ node1, node2 });
    }

    return node_q;
}

// 生成零位移边界值
vector<double> create_zero_u_values(vector<vector<int>> node_u)
{
    vector<double> values_u(node_u.size(), 0.0);

    return values_u;
}

// 生成均布边界载荷值
vector<vector<double>> create_uniform_q_values(
    vector<vector<int>> node_q,
    double qx,
    double qy)
{
    vector<vector<double>> values_q(
        node_q.size(),
        vector<double>(2, 0.0)
    );

    for (int i = 0; i < node_q.size(); i++)
    {
        values_q[i][0] = qx;
        values_q[i][1] = qy;
    }

    return values_q;
}

//六节点二次单元
vector<vector<int>> find_right_boundary_edges_T6(
    vector<vector<double>> P,
    vector<vector<int>> T6,
    double L)
{
    vector<vector<int>> node_q;
    double eps = 1.0e-8;

    for (int e = 0; e < T6.size(); e++)
    {
        int n1 = T6[e][0];
        int n2 = T6[e][1];
        int n3 = T6[e][2];
        int n4 = T6[e][3];
        int n5 = T6[e][4];
        int n6 = T6[e][5];

        vector<vector<int>> edges = {
            {n1, n4, n2},
            {n2, n5, n3},
            {n3, n6, n1}
        };

        for (int k = 0; k < edges.size(); k++)
        {
            int a = edges[k][0];
            int b = edges[k][1];
            int c = edges[k][2];

            double xa = P[a - 1][0];
            double xc = P[c - 1][0];

            if (abs(xa - L) < eps && abs(xc - L) < eps)
            {
                if (P[a - 1][1] < P[c - 1][1])
                {
                    node_q.push_back({ a, b, c });
                }
                else
                {
                    node_q.push_back({ c, b, a });
                }
            }
        }
    }

    sort(node_q.begin(), node_q.end(),
        [&](vector<int> edge1, vector<int> edge2)
        {
            double y1 = 0.5 * (P[edge1[0] - 1][1] + P[edge1[2] - 1][1]);
            double y2 = 0.5 * (P[edge2[0] - 1][1] + P[edge2[2] - 1][1]);

            return y1 < y2;
        });

    return node_q;
}

void apply_q_boundary_T6(
    vector<double>& F_global,
    vector<vector<double>> P,
    vector<vector<int>> node_q,
    vector<vector<double>> values_q,
    double thickness)
{
    for (int k = 0; k < node_q.size(); k++)
    {
        int node1 = node_q[k][0];
        int node2 = node_q[k][1];
        int node3 = node_q[k][2];

        double qx = values_q[k][0];
        double qy = values_q[k][1];

        double x1 = P[node1 - 1][0];
        double y1 = P[node1 - 1][1];

        double x3 = P[node3 - 1][0];
        double y3 = P[node3 - 1][1];

        double length = sqrt((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1));

        int dof1_x = 2 * (node1 - 1);
        int dof1_y = 2 * (node1 - 1) + 1;

        int dof2_x = 2 * (node2 - 1);
        int dof2_y = 2 * (node2 - 1) + 1;

        int dof3_x = 2 * (node3 - 1);
        int dof3_y = 2 * (node3 - 1) + 1;

        F_global[dof1_x] += qx * length * thickness / 6.0;
        F_global[dof1_y] += qy * length * thickness / 6.0;

        F_global[dof2_x] += qx * length * thickness * 4.0 / 6.0;
        F_global[dof2_y] += qy * length * thickness * 4.0 / 6.0;

        F_global[dof3_x] += qx * length * thickness / 6.0;
        F_global[dof3_y] += qy * length * thickness / 6.0;
    }
}