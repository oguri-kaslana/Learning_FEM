#include "fem_all.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace std;

void calculate_strain(const double B[3][8],
                      const double ue[8],
                      double strain[3])
{
    for (int i = 0; i < 3; ++i) {
        strain[i] = 0.0;
        for (int j = 0; j < 8; ++j) strain[i] += B[i][j] * ue[j];
    }
}

void calculate_stress(const MaterialParam& mat,
                      const double strain[3],
                      double stress[3])
{
    for (int i = 0; i < 3; ++i) {
        stress[i] = mat.D[i][0] * strain[0]
                  + mat.D[i][1] * strain[1]
                  + mat.D[i][2] * strain[2];
    }
}

double calculate_von_mises(double sigma_x,
                           double sigma_y,
                           double tau_xy)
{
    const double v = sigma_x * sigma_x - sigma_x * sigma_y + sigma_y * sigma_y
                   + 3.0 * tau_xy * tau_xy;
    return sqrt(max(0.0, v));
}

void calculate_element_center(const ElementInfo& elem,
                              const vector<NodeInfo>& nodes,
                              double& x_center,
                              double& y_center)
{
    x_center = 0.0;
    y_center = 0.0;
    for (int i = 0; i < 4; ++i) {
        const NodeInfo& nd = nodes[elem.node[i]];
        x_center += nd.x;
        y_center += nd.y;
    }
    x_center *= 0.25;
    y_center *= 0.25;
}

static ElementResult calculate_one_element_center(const ElementInfo& elem,
                                                  const vector<NodeInfo>& nodes,
                                                  const vector<double>& U,
                                                  const MaterialParam& mat,
                                                  bool use_fast)
{
    ElementResult er;
    er.elem_id = elem.id;
    calculate_element_center(elem, nodes, er.x_center, er.y_center);

    double ue[8];
    for (int i = 0; i < 8; ++i) ue[i] = U[elem.dof[i]];

    double B[3][8];
    double detJ = 0.0;
    create_B_center(elem, nodes, use_fast, B, detJ);

    double strain[3];
    double stress[3];
    calculate_strain(B, ue, strain);
    calculate_stress(mat, strain, stress);

    er.epsilon_x = strain[0];
    er.epsilon_y = strain[1];
    er.gamma_xy = strain[2];

    er.sigma_x = stress[0];
    er.sigma_y = stress[1];
    er.tau_xy = stress[2];
    er.sigma_vm = calculate_von_mises(er.sigma_x, er.sigma_y, er.tau_xy);

    return er;
}

void calculate_result_base(ModelData& model,
                           const MaterialParam& mat,
                           TimeRecord& timer,
                           const string& mode_name)
{
    timer.start_time("calculate_result_base");

    const vector<double>& U = model.U_ref.empty() ? model.U : model.U_ref;
    model.elem_result.assign(model.elements.size(), ElementResult());

    for (size_t e = 0; e < model.elements.size(); ++e) {
        model.elem_result[e] = calculate_one_element_center(model.elements[e], model.nodes, U, mat, false);
    }

    model.elem_result_base = model.elem_result;
    timer.stop_time("calculate_result_base");
    (void)mode_name;
}

void calculate_result_fast(ModelData& model,
                           const MaterialParam& mat,
                           TimeRecord& timer,
                           const string& mode_name)
{
    timer.start_time("calculate_result_fast");

    const vector<double>& U = model.U_ref.empty() ? model.U : model.U_ref;
    model.elem_result.assign(model.elements.size(), ElementResult());

#ifdef USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int e = 0; e < static_cast<int>(model.elements.size()); ++e) {
        model.elem_result[static_cast<size_t>(e)] =
            calculate_one_element_center(model.elements[static_cast<size_t>(e)], model.nodes, U, mat, true);
    }

    model.elem_result_fast = model.elem_result;
    timer.stop_time("calculate_result_fast");
    (void)mode_name;
}

void calculate_result_gauss_average(ModelData& model,
                                    const MaterialParam& mat,
                                    bool use_fast)
{
    const vector<double>& U = model.U_ref.empty() ? model.U : model.U_ref;
    model.elem_result.assign(model.elements.size(), ElementResult());

    if (model.gauss_data.empty()) {
        // 正常流程中 gauss_data 已由 create_Gauss 创建。这里保底使用中心点。
        for (size_t e = 0; e < model.elements.size(); ++e) {
            model.elem_result[e] = calculate_one_element_center(model.elements[e], model.nodes, U, mat, use_fast);
        }
        return;
    }

    for (size_t e = 0; e < model.elements.size(); ++e) {
        const ElementInfo& elem = model.elements[e];
        ElementResult er;
        er.elem_id = elem.id;
        calculate_element_center(elem, model.nodes, er.x_center, er.y_center);

        double ue[8];
        for (int i = 0; i < 8; ++i) ue[i] = U[elem.dof[i]];

        double strain_sum[3] = {0.0, 0.0, 0.0};
        double stress_sum[3] = {0.0, 0.0, 0.0};
        int n_ok = 0;

        for (size_t ig = 0; ig < model.gauss_data.size(); ++ig) {
            double B[3][8];
            double detJ = 0.0;
            create_B(elem, model.nodes, model.gauss_data[ig], use_fast, B, detJ);
            if (detJ <= 0.0) continue;

            double strain[3];
            double stress[3];
            calculate_strain(B, ue, strain);
            calculate_stress(mat, strain, stress);
            for (int k = 0; k < 3; ++k) {
                strain_sum[k] += strain[k];
                stress_sum[k] += stress[k];
            }
            n_ok++;
        }

        if (n_ok > 0) {
            er.epsilon_x = strain_sum[0] / n_ok;
            er.epsilon_y = strain_sum[1] / n_ok;
            er.gamma_xy = strain_sum[2] / n_ok;
            er.sigma_x = stress_sum[0] / n_ok;
            er.sigma_y = stress_sum[1] / n_ok;
            er.tau_xy = stress_sum[2] / n_ok;
            er.sigma_vm = calculate_von_mises(er.sigma_x, er.sigma_y, er.tau_xy);
        }

        model.elem_result[e] = er;
    }
}

void calculate_summary(ModelData& model,
                       const GeometryParam& geom,
                       const LoadParam& load)
{
    CaseSummary s;
    s.ratio = geom.ratio;
    s.a = geom.a;
    s.b = geom.b;

    const vector<double>& U = model.U_ref.empty() ? model.U : model.U_ref;
    for (size_t i = 0; i + 1 < U.size(); i += 2) {
        const double ux = U[i];
        const double uy = U[i + 1];
        const double umag = sqrt(ux * ux + uy * uy);
        s.max_ux = max(s.max_ux, fabs(ux));
        s.max_uy = max(s.max_uy, fabs(uy));
        s.max_umag = max(s.max_umag, umag);
    }

    s.max_sigma_x = 0.0;
    s.max_sigma_y = 0.0;
    s.max_tau_xy = 0.0;
    s.max_sigma_vm = 0.0;

    for (size_t i = 0; i < model.elem_result.size(); ++i) {
        const ElementResult& er = model.elem_result[i];
        if (fabs(er.sigma_x) > fabs(s.max_sigma_x)) s.max_sigma_x = er.sigma_x;
        if (fabs(er.sigma_y) > fabs(s.max_sigma_y)) s.max_sigma_y = er.sigma_y;
        if (fabs(er.tau_xy) > fabs(s.max_tau_xy)) s.max_tau_xy = er.tau_xy;

        if (er.sigma_vm > s.max_sigma_vm) {
            s.max_sigma_vm = er.sigma_vm;
            s.max_sigma_vm_elem = er.elem_id;
            s.max_sigma_vm_x = er.x_center;
            s.max_sigma_vm_y = er.y_center;
        }
    }

    const double q = max(1.0e-30, fabs(load.q));
    s.Kt_x = fabs(s.max_sigma_x) / q;
    s.Kt_vm = s.max_sigma_vm / q;

    model.summary = s;
}

void calculate_result(ModelData& model,
                      const MaterialParam& mat,
                      const LoadParam& load,
                      const SpeedParam& speed,
                      TimeRecord& timer,
                      const string& mode_name)
{
    if (speed.stress_mode == CENTER_STRESS) {
        if (speed.use_fast_post) calculate_result_fast(model, mat, timer, mode_name);
        else calculate_result_base(model, mat, timer, mode_name);
    } else if (speed.stress_mode == GAUSS_AVERAGE_STRESS) {
        timer.start_time(speed.use_fast_post ? "calculate_result_fast" : "calculate_result_base");
        calculate_result_gauss_average(model, mat, speed.use_fast_post);
        timer.stop_time(speed.use_fast_post ? "calculate_result_fast" : "calculate_result_base");
    } else {
        calculate_result_base(model, mat, timer, mode_name);
        vector<ElementResult> center_result = model.elem_result;
        timer.start_time("calculate_result_fast");
        calculate_result_gauss_average(model, mat, true);
        timer.stop_time("calculate_result_fast");
        model.elem_result_base = center_result;
        model.elem_result_fast = model.elem_result;
    }

    // 工况汇总由 main 中的 calculate_summary(model, geom, load) 单独调用，避免几何参数传递混乱。
    (void)load;
    (void)mode_name;
}
