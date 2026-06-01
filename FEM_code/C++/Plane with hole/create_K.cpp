#include "fem_all.h"

#include <iostream>
#include <cmath>
#include <algorithm>

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace std;

void clear_local_K(double Ke[8][8])
{
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            Ke[i][j] = 0.0;
        }
    }
}

static void multiply_D_B(const MaterialParam& mat,
                         const double B[3][8],
                         double DB[3][8])
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 8; ++j) {
            DB[i][j] = mat.D[i][0] * B[0][j]
                     + mat.D[i][1] * B[1][j]
                     + mat.D[i][2] * B[2][j];
        }
    }
}

void create_K_local(const ElementInfo& elem,
                    const vector<NodeInfo>& nodes,
                    const MaterialParam& mat,
                    const vector<GaussData>& gauss_data,
                    bool use_fast,
                    double Ke[8][8])
{
    clear_local_K(Ke);

    if (gauss_data.empty()) {
        cerr << "[create_K warning] gauss_data is empty." << endl;
        return;
    }

    for (size_t ig = 0; ig < gauss_data.size(); ++ig) {
        double B[3][8];
        double DB[3][8];
        double detJ = 0.0;

        create_B(elem, nodes, gauss_data[ig], use_fast, B, detJ);

        if (detJ <= 0.0) {
            // detJ 非正时仍继续计算会污染结果，因此这里跳过该高斯点。
            continue;
        }

        multiply_D_B(mat, B, DB);

        const double coef = mat.t * detJ * gauss_data[ig].weight;
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                double v = 0.0;
                for (int k = 0; k < 3; ++k) {
                    v += B[k][i] * DB[k][j];
                }
                Ke[i][j] += v * coef;
            }
        }
    }
}

void create_K_dense_base(ModelData& model,
                         const MaterialParam& mat,
                         TimeRecord& timer,
                         const string& mode_name)
{
    const string step = "create_K_dense_base";
    timer.start_time(step);

    model.n_node = static_cast<int>(model.nodes.size());
    model.n_elem = static_cast<int>(model.elements.size());
    model.n_dof = 2 * model.n_node;
    model.K_dense.assign(static_cast<size_t>(model.n_dof) * static_cast<size_t>(model.n_dof), 0.0);

    for (size_t e = 0; e < model.elements.size(); ++e) {
        const ElementInfo& elem = model.elements[e];
        double Ke[8][8];
        create_K_local(elem, model.nodes, mat, model.gauss_data, false, Ke);

        for (int i = 0; i < 8; ++i) {
            const int row = elem.dof[i];
            for (int j = 0; j < 8; ++j) {
                const int col = elem.dof[j];
                model.K_dense[static_cast<size_t>(row) * model.n_dof + col] += Ke[i][j];
            }
        }
    }

    timer.stop_time(step);
    (void)mode_name;
}

void create_K_dense_fast(ModelData& model,
                         const MaterialParam& mat,
                         TimeRecord& timer,
                         const string& mode_name)
{
    const string step = "create_K_dense_fast";
    timer.start_time(step);

    model.n_node = static_cast<int>(model.nodes.size());
    model.n_elem = static_cast<int>(model.elements.size());
    model.n_dof = 2 * model.n_node;
    model.K_dense.assign(static_cast<size_t>(model.n_dof) * static_cast<size_t>(model.n_dof), 0.0);

    // dense 全局矩阵装配存在写冲突，默认串行保证结果稳定。
    for (size_t e = 0; e < model.elements.size(); ++e) {
        const ElementInfo& elem = model.elements[e];
        double Ke[8][8];
        create_K_local(elem, model.nodes, mat, model.gauss_data, true, Ke);

        const int* edof = elem.dof;
        for (int i = 0; i < 8; ++i) {
            const int row = edof[i];
            const size_t base = static_cast<size_t>(row) * static_cast<size_t>(model.n_dof);
            for (int j = 0; j < 8; ++j) {
                model.K_dense[base + edof[j]] += Ke[i][j];
            }
        }
    }

    timer.stop_time(step);
    (void)mode_name;
}

void create_K_sparse_coo(ModelData& model,
                         const MaterialParam& mat,
                         TimeRecord& timer,
                         const string& mode_name)
{
    const string step = "create_K_sparse_coo";
    timer.start_time(step);

    model.n_node = static_cast<int>(model.nodes.size());
    model.n_elem = static_cast<int>(model.elements.size());
    model.n_dof = 2 * model.n_node;

    model.K_coo.resize_matrix(model.n_dof, model.n_dof);
    model.K_coo.reserve_data(static_cast<size_t>(model.elements.size()) * 64u);

    for (size_t e = 0; e < model.elements.size(); ++e) {
        const ElementInfo& elem = model.elements[e];
        double Ke[8][8];
        create_K_local(elem, model.nodes, mat, model.gauss_data, true, Ke);

        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                if (fabs(Ke[i][j]) > 0.0) {
                    model.K_coo.add_value(elem.dof[i], elem.dof[j], Ke[i][j]);
                }
            }
        }
    }

    timer.stop_time(step);

    timer.start_time("convert_coo_to_csr");
    convert_coo_csr(model.K_coo, model.K_csr);
    model.sparse_nnz = model.K_csr.get_nnz();
    model.sparse_rate = calculate_sparse_rate(model.K_csr);
    timer.stop_time("convert_coo_to_csr");

    (void)mode_name;
}

void create_K_sparse_parallel(ModelData& model,
                              const MaterialParam& mat,
                              const SpeedParam& speed,
                              TimeRecord& timer,
                              const string& mode_name)
{
    const string step = "create_K_sparse_parallel";
    timer.start_time(step);

    // 默认实现为安全串行 COO；如果启用 OpenMP，使用线程局部 COO 后合并。
    model.n_node = static_cast<int>(model.nodes.size());
    model.n_elem = static_cast<int>(model.elements.size());
    model.n_dof = 2 * model.n_node;
    model.K_coo.resize_matrix(model.n_dof, model.n_dof);

#ifdef USE_OPENMP
    if (speed.use_parallel_loop) {
        int n_thread = omp_get_max_threads();
        vector<SparseCOO> local_coo(n_thread);
        for (int k = 0; k < n_thread; ++k) {
            local_coo[k].resize_matrix(model.n_dof, model.n_dof);
            local_coo[k].reserve_data(static_cast<size_t>(model.elements.size()) * 64u / n_thread + 64u);
        }

#pragma omp parallel for schedule(static)
        for (int e = 0; e < static_cast<int>(model.elements.size()); ++e) {
            const int tid = omp_get_thread_num();
            const ElementInfo& elem = model.elements[static_cast<size_t>(e)];
            double Ke[8][8];
            create_K_local(elem, model.nodes, mat, model.gauss_data, true, Ke);

            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < 8; ++j) {
                    if (fabs(Ke[i][j]) > 0.0) {
                        local_coo[tid].add_value(elem.dof[i], elem.dof[j], Ke[i][j]);
                    }
                }
            }
        }

        size_t total_raw = 0;
        for (int k = 0; k < n_thread; ++k) total_raw += local_coo[k].data.size();
        model.K_coo.reserve_data(total_raw);
        for (int k = 0; k < n_thread; ++k) {
            model.K_coo.data.insert(model.K_coo.data.end(),
                                    local_coo[k].data.begin(),
                                    local_coo[k].data.end());
        }
    } else
#endif
    {
        model.K_coo.reserve_data(static_cast<size_t>(model.elements.size()) * 64u);
        for (size_t e = 0; e < model.elements.size(); ++e) {
            const ElementInfo& elem = model.elements[e];
            double Ke[8][8];
            create_K_local(elem, model.nodes, mat, model.gauss_data, true, Ke);
            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < 8; ++j) {
                    if (fabs(Ke[i][j]) > 0.0) model.K_coo.add_value(elem.dof[i], elem.dof[j], Ke[i][j]);
                }
            }
        }
    }

    timer.stop_time(step);

    timer.start_time("convert_coo_to_csr");
    convert_coo_csr(model.K_coo, model.K_csr);
    model.sparse_nnz = model.K_csr.get_nnz();
    model.sparse_rate = calculate_sparse_rate(model.K_csr);
    timer.stop_time("convert_coo_to_csr");

    (void)mode_name;
}

void create_K_global(ModelData& model,
                     const MaterialParam& mat,
                     const SpeedParam& speed,
                     TimeRecord& timer,
                     const string& mode_name)
{
    // 确保高斯点存在。
    if (model.gauss_data.empty()) {
        create_Gauss(model.gauss_data, timer, mode_name);
    }

    if (speed.matrix_mode == DENSE_MATRIX || speed.matrix_mode == BOTH_MATRIX) {
        if (speed.assemble_mode == ASSEMBLE_DENSE_BASE ||
            speed.assemble_mode == ASSEMBLE_ALL ||
            !speed.use_fast_assemble) {
            create_K_dense_base(model, mat, timer, mode_name);
        }

        if (speed.assemble_mode == ASSEMBLE_DENSE_FAST ||
            speed.assemble_mode == ASSEMBLE_ALL ||
            speed.use_fast_assemble) {
            create_K_dense_fast(model, mat, timer, mode_name);
        }
    }

    if (speed.matrix_mode == SPARSE_MATRIX || speed.matrix_mode == BOTH_MATRIX || speed.use_sparse_matrix) {
        if (speed.use_parallel_loop) {
            create_K_sparse_parallel(model, mat, speed, timer, mode_name);
        } else {
            create_K_sparse_coo(model, mat, timer, mode_name);
        }
    }
}
