#include "fem_all.h"

#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

static void add_fixed_dof(ModelData& model, int dof)
{
    if (dof < 0 || dof >= model.n_dof) return;
    if (model.is_fixed.empty()) model.is_fixed.assign(model.n_dof, 0);
    if (!model.is_fixed[dof]) {
        model.is_fixed[dof] = 1;
        model.fixed_dof.push_back(dof);
    }
}

void apply_fixed_boundary(ModelData& model)
{
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        if (model.nodes[i].flag_left == 1) {
            add_fixed_dof(model, 2 * model.nodes[i].id);
            add_fixed_dof(model, 2 * model.nodes[i].id + 1);
        }
    }
}

void apply_symmetry_boundary(ModelData& model)
{
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        if (model.nodes[i].flag_sym == 1) {
            // 水平中线对称：只约束 v = 0，不约束 u。
            add_fixed_dof(model, 2 * model.nodes[i].id + 1);
        }
    }
}

void create_dof_set(ModelData& model,
                    TimeRecord& timer,
                    const string& mode_name)
{
    timer.start_time("create_free_dof");

    model.n_node = static_cast<int>(model.nodes.size());
    model.n_elem = static_cast<int>(model.elements.size());
    model.n_dof = 2 * model.n_node;

    if (static_cast<int>(model.is_fixed.size()) != model.n_dof) {
        model.is_fixed.assign(model.n_dof, 0);
    }

    model.fixed_dof.clear();
    model.free_dof.clear();

    apply_fixed_boundary(model);
    apply_symmetry_boundary(model);

    sort(model.fixed_dof.begin(), model.fixed_dof.end());
    model.fixed_dof.erase(unique(model.fixed_dof.begin(), model.fixed_dof.end()), model.fixed_dof.end());

    model.free_dof.reserve(model.n_dof - static_cast<int>(model.fixed_dof.size()));
    for (int i = 0; i < model.n_dof; ++i) {
        if (!model.is_fixed[i]) model.free_dof.push_back(i);
    }

    model.n_left_fixed_dof = 0;
    model.n_sym_fixed_dof = 0;
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        if (model.nodes[i].flag_left) model.n_left_fixed_dof += 2;
        if (model.nodes[i].flag_sym) model.n_sym_fixed_dof += 1;
    }

    timer.stop_time("create_free_dof");
    (void)mode_name;
}

void apply_traction_load(ModelData& model,
                         const GeometryParam& geom,
                         LoadParam& load,
                         TimeRecord& timer,
                         const string& mode_name)
{
    timer.start_time("apply_load");

    if (model.n_dof <= 0) {
        model.n_node = static_cast<int>(model.nodes.size());
        model.n_dof = 2 * model.n_node;
    }
    if (static_cast<int>(model.F.size()) != model.n_dof) {
        model.F.assign(model.n_dof, 0.0);
    }

    load.total_force = 0.0;
    model.n_right_load_edge = 0;

    for (size_t k = 0; k < model.edges.size(); ++k) {
        const BoundaryEdge& edge = model.edges[k];
        if (edge.marker != 2) continue;

        const NodeInfo& n1 = model.nodes[edge.n1];
        const NodeInfo& n2 = model.nodes[edge.n2];

        // 只对 x = L 的右边界施加载荷，避免误标边界。
        if (fabs(n1.x - geom.L) > 1.0e-8 || fabs(n2.x - geom.L) > 1.0e-8) continue;

        const double dx = n2.x - n1.x;
        const double dy = n2.y - n1.y;
        const double le = sqrt(dx * dx + dy * dy);
        const double f_node = load.q * le * 0.5; // t 在 MaterialParam 中，这里默认 t=1；统一接口下 check 使用 mat.t。

        if (load.direction == 0) {
            model.F[2 * n1.id] += f_node;
            model.F[2 * n2.id] += f_node;
        } else {
            model.F[2 * n1.id + 1] += f_node;
            model.F[2 * n2.id + 1] += f_node;
        }

        load.total_force += 2.0 * f_node;
        model.n_right_load_edge += 1;
    }

    model.total_load_actual = load.total_force;
    model.total_load_theory = load.q * (geom.H / 2.0); // t=1 时的理论半模型外力

    timer.stop_time("apply_load");
    (void)mode_name;
}

void apply_load(ModelData& model,
                const GeometryParam& geom,
                LoadParam& load,
                TimeRecord& timer,
                const string& mode_name)
{
    apply_traction_load(model, geom, load, timer, mode_name);
}

static void apply_boundary_dense_core(ModelData& model)
{
    if (model.K_dense.empty() || model.n_dof <= 0) return;

    const int n = model.n_dof;

    // 这里已知位移均为 0，因此无需执行 F -= K(:,id)*u_known。
    for (size_t k = 0; k < model.fixed_dof.size(); ++k) {
        const int id = model.fixed_dof[k];
        if (id < 0 || id >= n) continue;

        for (int j = 0; j < n; ++j) {
            model.K_dense[static_cast<size_t>(id) * n + j] = 0.0;
        }
        for (int i = 0; i < n; ++i) {
            model.K_dense[static_cast<size_t>(i) * n + id] = 0.0;
        }
        model.K_dense[static_cast<size_t>(id) * n + id] = 1.0;
        model.F[id] = 0.0;
    }
}

void apply_boundary_dense_base(ModelData& model,
                               TimeRecord& timer,
                               const string& mode_name)
{
    timer.start_time("apply_boundary_base");
    if (model.free_dof.empty() && model.fixed_dof.empty()) {
        create_dof_set(model, timer, mode_name);
    }
    apply_boundary_dense_core(model);
    timer.stop_time("apply_boundary_base");
}

void apply_boundary_dense_fast(ModelData& model,
                               TimeRecord& timer,
                               const string& mode_name)
{
    timer.start_time("apply_boundary_fast");
    if (model.free_dof.empty() && model.fixed_dof.empty()) {
        create_dof_set(model, timer, mode_name);
    }
    apply_boundary_dense_core(model);
    timer.stop_time("apply_boundary_fast");
}

void apply_boundary_sparse(ModelData& model,
                           TimeRecord& timer,
                           const string& mode_name)
{
    timer.start_time("apply_boundary_sparse");

    if (model.free_dof.empty() && model.fixed_dof.empty()) {
        create_dof_set(model, timer, mode_name);
    }

    if (model.K_csr.n_row > 0) {
        extract_free_csr(model.K_csr, model.free_dof, model.Kff_csr);
        extract_free_vector(model.F, model.free_dof, model.Ff);
        model.Uf.assign(model.free_dof.size(), 0.0);
    }

    timer.stop_time("apply_boundary_sparse");
}

void apply_boundary(ModelData& model,
                    const SpeedParam& speed,
                    TimeRecord& timer,
                    const string& mode_name)
{
    create_dof_set(model, timer, mode_name);

    if (speed.matrix_mode == DENSE_MATRIX || speed.matrix_mode == BOTH_MATRIX) {
        if (speed.use_fast_boundary) apply_boundary_dense_fast(model, timer, mode_name);
        else apply_boundary_dense_base(model, timer, mode_name);
    }

    if (speed.matrix_mode == SPARSE_MATRIX || speed.matrix_mode == BOTH_MATRIX || speed.use_sparse_matrix) {
        apply_boundary_sparse(model, timer, mode_name);
    }
}

bool check_total_load(ModelData& model,
                      const GeometryParam& geom,
                      const MaterialParam& mat,
                      const LoadParam& load,
                      double eps)
{
    model.total_load_theory = load.q * (geom.H / 2.0) * mat.t;

    double actual = 0.0;
    for (int i = 0; i < model.n_dof; i += 2) {
        actual += model.F[i];
    }
    model.total_load_actual = actual;

    const double den = max(1.0, fabs(model.total_load_theory));
    const double rel = fabs(actual - model.total_load_theory) / den;
    if (rel > eps) {
        model.add_warning("total load check failed: actual != theory");
        return false;
    }
    return true;
}
