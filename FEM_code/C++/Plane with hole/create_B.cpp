#include "fem_all.h"

#include <iostream>
#include <cmath>

using namespace std;

void clear_B(double B[3][8])
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 8; ++j) {
            B[i][j] = 0.0;
        }
    }
}

void calculate_dN_xy(const GaussData& gp,
                     const double invJ[2][2],
                     double dN_dx[4],
                     double dN_dy[4])
{
    // 由于 J 定义为：
    // [dx/dxi   dy/dxi]
    // [dx/deta  dy/deta]
    // 因此 [dN/dx dN/dy] = [dN/dxi dN/deta] * inv(J)
    for (int i = 0; i < 4; ++i) {
        dN_dx[i] = gp.dN_dxi[i] * invJ[0][0] + gp.dN_deta[i] * invJ[1][0];
        dN_dy[i] = gp.dN_dxi[i] * invJ[0][1] + gp.dN_deta[i] * invJ[1][1];
    }
}

static void fill_B_matrix(const double dN_dx[4],
                          const double dN_dy[4],
                          double B[3][8])
{
    clear_B(B);

    for (int i = 0; i < 4; ++i) {
        const int c = 2 * i;

        B[0][c]     = dN_dx[i];
        B[1][c + 1] = dN_dy[i];

        B[2][c]     = dN_dy[i];
        B[2][c + 1] = dN_dx[i];
    }
}

void create_B_base(const ElementInfo& elem,
                   const vector<NodeInfo>& nodes,
                   const GaussData& gp,
                   double B[3][8],
                   double& detJ)
{
    double J[2][2];
    double invJ[2][2];
    bool ok = create_Jacobi_base(elem, nodes, gp, J, invJ, detJ);

    if (!ok || detJ <= 0.0) {
        cerr << "[create_B warning] detJ <= 0 or singular Jacobi in element "
             << elem.id << ", detJ = " << detJ << endl;
    }

    double dN_dx[4];
    double dN_dy[4];
    calculate_dN_xy(gp, invJ, dN_dx, dN_dy);
    fill_B_matrix(dN_dx, dN_dy, B);
}

void create_B_fast(const ElementInfo& elem,
                   const vector<NodeInfo>& nodes,
                   const GaussData& gp,
                   double B[3][8],
                   double& detJ)
{
    double J[2][2];
    double invJ[2][2];
    bool ok = create_Jacobi_fast(elem, nodes, gp, J, invJ, detJ);

    if (!ok || detJ <= 0.0) {
        cerr << "[create_B warning] detJ <= 0 or singular Jacobi in element "
             << elem.id << ", detJ = " << detJ << endl;
    }

    double dN_dx[4];
    double dN_dy[4];
    calculate_dN_xy(gp, invJ, dN_dx, dN_dy);
    fill_B_matrix(dN_dx, dN_dy, B);
}

void create_B(const ElementInfo& elem,
              const vector<NodeInfo>& nodes,
              const GaussData& gp,
              bool use_fast,
              double B[3][8],
              double& detJ)
{
    if (use_fast) {
        create_B_fast(elem, nodes, gp, B, detJ);
    } else {
        create_B_base(elem, nodes, gp, B, detJ);
    }
}

void create_B_center(const ElementInfo& elem,
                     const vector<NodeInfo>& nodes,
                     bool use_fast,
                     double B[3][8],
                     double& detJ)
{
    GaussData gp;
    gp.xi = 0.0;
    gp.eta = 0.0;
    gp.weight = 4.0;
    create_shape_Q4(0.0, 0.0, gp.N, gp.dN_dxi, gp.dN_deta);
    create_B(elem, nodes, gp, use_fast, B, detJ);
}
