#include "fem_all.h"

#include <cmath>
#include <iostream>
#include <limits>

using namespace std;

void create_Gauss(vector<GaussData>& gauss_data,
                  TimeRecord& timer,
                  const string& mode_name)
{
    const string step = "create_Gauss";
    timer.start_time(step);

    gauss_data.clear();
    gauss_data.reserve(4);

    const double g = 1.0 / sqrt(3.0);
    const double xi_list[2] = {-g, g};
    const double eta_list[2] = {-g, g};

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            GaussData gp;
            gp.xi = xi_list[j];
            gp.eta = eta_list[i];
            gp.weight = 1.0;
            create_shape_Q4(gp.xi, gp.eta, gp.N, gp.dN_dxi, gp.dN_deta);
            gauss_data.push_back(gp);
        }
    }

    timer.stop_time(step);
    (void)mode_name;
}

void create_shape_Q4(double xi,
                     double eta,
                     double N[4],
                     double dN_dxi[4],
                     double dN_deta[4])
{
    // Q4 节点顺序：左下、右下、右上、左上，对应自然坐标逆时针。
    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);

    dN_dxi[0] = -0.25 * (1.0 - eta);
    dN_dxi[1] =  0.25 * (1.0 - eta);
    dN_dxi[2] =  0.25 * (1.0 + eta);
    dN_dxi[3] = -0.25 * (1.0 + eta);

    dN_deta[0] = -0.25 * (1.0 - xi);
    dN_deta[1] = -0.25 * (1.0 + xi);
    dN_deta[2] =  0.25 * (1.0 + xi);
    dN_deta[3] =  0.25 * (1.0 - xi);
}

bool inverse_matrix_2x2(const double J[2][2],
                        double invJ[2][2],
                        double& detJ)
{
    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

    if (fabs(detJ) < 1.0e-14) {
        invJ[0][0] = invJ[0][1] = invJ[1][0] = invJ[1][1] = 0.0;
        return false;
    }

    invJ[0][0] =  J[1][1] / detJ;
    invJ[0][1] = -J[0][1] / detJ;
    invJ[1][0] = -J[1][0] / detJ;
    invJ[1][1] =  J[0][0] / detJ;

    return true;
}

bool create_Jacobi_base(const ElementInfo& elem,
                        const vector<NodeInfo>& nodes,
                        const GaussData& gp,
                        double J[2][2],
                        double invJ[2][2],
                        double& detJ)
{
    J[0][0] = 0.0;  // dx/dxi
    J[0][1] = 0.0;  // dy/dxi
    J[1][0] = 0.0;  // dx/deta
    J[1][1] = 0.0;  // dy/deta

    for (int i = 0; i < 4; ++i) {
        const NodeInfo& nd = nodes[elem.node[i]];
        J[0][0] += gp.dN_dxi[i]  * nd.x;
        J[0][1] += gp.dN_dxi[i]  * nd.y;
        J[1][0] += gp.dN_deta[i] * nd.x;
        J[1][1] += gp.dN_deta[i] * nd.y;
    }

    return inverse_matrix_2x2(J, invJ, detJ);
}

bool create_Jacobi_fast(const ElementInfo& elem,
                        const vector<NodeInfo>& nodes,
                        const GaussData& gp,
                        double J[2][2],
                        double invJ[2][2],
                        double& detJ)
{
    // 加速版仍保持清晰写法，区别在于后续 create_K 中会复用固定数组并避免 vector 临时对象。
    const NodeInfo& n1 = nodes[elem.node[0]];
    const NodeInfo& n2 = nodes[elem.node[1]];
    const NodeInfo& n3 = nodes[elem.node[2]];
    const NodeInfo& n4 = nodes[elem.node[3]];

    const double x[4] = {n1.x, n2.x, n3.x, n4.x};
    const double y[4] = {n1.y, n2.y, n3.y, n4.y};

    J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;

    for (int i = 0; i < 4; ++i) {
        J[0][0] += gp.dN_dxi[i]  * x[i];
        J[0][1] += gp.dN_dxi[i]  * y[i];
        J[1][0] += gp.dN_deta[i] * x[i];
        J[1][1] += gp.dN_deta[i] * y[i];
    }

    return inverse_matrix_2x2(J, invJ, detJ);
}

bool create_Jacobi(const ElementInfo& elem,
                   const vector<NodeInfo>& nodes,
                   const GaussData& gp,
                   bool use_fast,
                   double J[2][2],
                   double invJ[2][2],
                   double& detJ)
{
    if (use_fast) {
        return create_Jacobi_fast(elem, nodes, gp, J, invJ, detJ);
    }
    return create_Jacobi_base(elem, nodes, gp, J, invJ, detJ);
}

double calculate_detJ(const ElementInfo& elem,
                      const vector<NodeInfo>& nodes,
                      double xi,
                      double eta)
{
    GaussData gp;
    gp.xi = xi;
    gp.eta = eta;
    gp.weight = 1.0;
    create_shape_Q4(xi, eta, gp.N, gp.dN_dxi, gp.dN_deta);

    double J[2][2];
    double invJ[2][2];
    double detJ = 0.0;
    create_Jacobi_fast(elem, nodes, gp, J, invJ, detJ);
    return detJ;
}
