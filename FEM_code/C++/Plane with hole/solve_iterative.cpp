#include "fem_all.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;

static SolverInfo make_iter_info(const string& method,
                                 const string& matrix_type,
                                 double omega,
                                 int n_dof,
                                 int n_solve_dof)
{
    SolverInfo info;
    info.method = method;
    info.matrix_type = matrix_type;
    info.omega = omega;
    info.converged = false;
    info.iterations = 0;
    info.time_seconds = 0.0;
    info.final_residual = numeric_limits<double>::max();
    info.max_displacement = 0.0;
    info.n_dof = n_dof;
    info.n_solve_dof = n_solve_dof;
    return info;
}

static void prepare_dense_free(const ModelData& model,
                               vector<double>& A,
                               vector<double>& b,
                               int& n)
{
    extract_free_dof_dense(model.K_dense, model.F, model.free_dof, model.n_dof, A, b);
    n = static_cast<int>(model.free_dof.size());
}

static SparseCSR prepare_sparse_free(ModelData& model,
                                     vector<double>& b)
{
    if (model.Kff_csr.n_row > 0) {
        b = model.Ff;
        return model.Kff_csr;
    }
    SparseCSR Kff;
    extract_free_dof_sparse(model.K_csr, model.F, model.free_dof, Kff, b);
    return Kff;
}

static void finish_iter(ModelData& model,
                        const vector<double>& x_free,
                        SolverInfo& info,
                        double residual)
{
    recover_full_vector(x_free, model.free_dof, model.n_dof, model.U);
    info.final_residual = residual;
    info.max_displacement = calculate_max_displacement(model.U);
    model.iterative_solver_info.push_back(info);
    model.all_solver_info.push_back(info);
}

static void record_iter_history(IterHistory* hist,
                                int iter,
                                double residual,
                                const chrono::high_resolution_clock::time_point& t0)
{
    if (hist == NULL) return;

    // 前 200 步全部记录，之后每 10 步记录一次；最终收敛步由外层调用补充。
    // 这样既能画出曲线，又避免 Jacobi / GS / SOR 产生过大的 csv 文件。
    if (iter <= 200 || iter % 10 == 0) {
        const auto now = chrono::high_resolution_clock::now();
        const double t = chrono::duration<double>(now - t0).count();
        hist->iteration.push_back(iter);
        hist->residual.push_back(residual);
        hist->time_seconds.push_back(t);
    }
}

static void record_iter_final(IterHistory* hist,
                              int iter,
                              double residual,
                              const chrono::high_resolution_clock::time_point& t0)
{
    if (hist == NULL) return;
    if (!hist->iteration.empty() && hist->iteration.back() == iter) return;
    const auto now = chrono::high_resolution_clock::now();
    const double t = chrono::duration<double>(now - t0).count();
    hist->iteration.push_back(iter);
    hist->residual.push_back(residual);
    hist->time_seconds.push_back(t);
}

static void init_iter_history(IterHistory& hist, const SolverInfo& info)
{
    hist.method = info.method;
    hist.matrix_type = info.matrix_type;
    hist.omega = info.omega;
    hist.iteration.clear();
    hist.residual.clear();
    hist.time_seconds.clear();
}

static void push_iter_history(ModelData& model, const IterHistory& hist)
{
    if (!hist.iteration.empty()) {
        model.iterative_history.push_back(hist);
    }
}

static bool iterate_jacobi_dense_core(const vector<double>& A,
                                      const vector<double>& b,
                                      int n,
                                      double tol,
                                      int max_iter,
                                      vector<double>& x,
                                      int& iter,
                                      double& res,
                                      IterHistory* hist)
{
    x.assign(n, 0.0);
    vector<double> x_old(n, 0.0);
    const auto t0 = chrono::high_resolution_clock::now();

    for (iter = 1; iter <= max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double diag = A[static_cast<size_t>(i) * n + i];
            if (fabs(diag) < 1.0e-30) return false;
            double s = b[i];
            for (int j = 0; j < n; ++j) {
                if (j != i) s -= A[static_cast<size_t>(i) * n + j] * x_old[j];
            }
            x[i] = s / diag;
        }
        res = calculate_residual_dense(A, x, b, n);
        record_iter_history(hist, iter, res, t0);
        if (res < tol) { record_iter_final(hist, iter, res, t0); return true; }
        x_old.swap(x);
    }
    record_iter_final(hist, iter > max_iter ? max_iter : iter, res, t0);
    return false;
}

static bool iterate_gs_dense_core(const vector<double>& A,
                                  const vector<double>& b,
                                  int n,
                                  double tol,
                                  int max_iter,
                                  double omega,
                                  vector<double>& x,
                                  int& iter,
                                  double& res,
                                  IterHistory* hist)
{
    x.assign(n, 0.0);
    const auto t0 = chrono::high_resolution_clock::now();
    for (iter = 1; iter <= max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double diag = A[static_cast<size_t>(i) * n + i];
            if (fabs(diag) < 1.0e-30) return false;
            double s = b[i];
            for (int j = 0; j < n; ++j) {
                if (j != i) s -= A[static_cast<size_t>(i) * n + j] * x[j];
            }
            const double x_gs = s / diag;
            x[i] = (1.0 - omega) * x[i] + omega * x_gs;
        }
        res = calculate_residual_dense(A, x, b, n);
        record_iter_history(hist, iter, res, t0);
        if (res < tol) { record_iter_final(hist, iter, res, t0); return true; }
        if (!isfinite(res) || res > 1.0e80) return false;
    }
    record_iter_final(hist, iter > max_iter ? max_iter : iter, res, t0);
    return false;
}

static bool iterate_cg_dense_core(const vector<double>& A,
                                  const vector<double>& b,
                                  int n,
                                  double tol,
                                  int max_iter,
                                  bool use_pcg,
                                  vector<double>& x,
                                  int& iter,
                                  double& res,
                                  IterHistory* hist)
{
    x.assign(n, 0.0);
    const auto t0 = chrono::high_resolution_clock::now();
    vector<double> r = b;
    vector<double> z(n, 0.0), p(n, 0.0), Ap(n, 0.0), diag(n, 1.0);

    if (use_pcg) {
        extract_diag_dense(A, n, diag);
        for (int i = 0; i < n; ++i) {
            if (fabs(diag[i]) < 1.0e-30) diag[i] = 1.0;
            z[i] = r[i] / diag[i];
        }
        p = z;
    } else {
        p = r;
        z = r;
    }

    double rz_old = dot_vector(r, z);
    const double nb = max(1.0e-30, norm_vector(b));
    res = norm_vector(r) / nb;
    record_iter_history(hist, 0, res, t0);
    if (res < tol) { iter = 0; record_iter_final(hist, iter, res, t0); return true; }

    for (iter = 1; iter <= max_iter; ++iter) {
        multiply_matrix_vector_dense(A, p, Ap, n);
        const double denom = dot_vector(p, Ap);
        if (fabs(denom) < 1.0e-30) return false;
        const double alpha = rz_old / denom;

        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        res = norm_vector(r) / nb;
        record_iter_history(hist, iter, res, t0);
        if (res < tol) { record_iter_final(hist, iter, res, t0); return true; }

        if (use_pcg) {
            for (int i = 0; i < n; ++i) z[i] = r[i] / diag[i];
        } else {
            z = r;
        }

        const double rz_new = dot_vector(r, z);
        const double beta = rz_new / max(1.0e-300, rz_old);
        for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        rz_old = rz_new;
    }
    record_iter_final(hist, iter > max_iter ? max_iter : iter, res, t0);
    return false;
}

static bool iterate_jacobi_sparse_core(const SparseCSR& A,
                                       const vector<double>& b,
                                       double tol,
                                       int max_iter,
                                       vector<double>& x,
                                       int& iter,
                                       double& res,
                                       IterHistory* hist)
{
    const int n = A.n_row;
    x.assign(n, 0.0);
    vector<double> x_old(n, 0.0);
    const auto t0 = chrono::high_resolution_clock::now();
    vector<double> diag(n, 0.0);
    extract_diag_sparse(A, diag);

    for (iter = 1; iter <= max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            if (fabs(diag[i]) < 1.0e-30) return false;
            double s = b[i];
            for (int k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
                int j = A.col_id[k];
                if (j != i) s -= A.val[k] * x_old[j];
            }
            x[i] = s / diag[i];
        }
        res = calculate_residual_sparse(A, x, b);
        record_iter_history(hist, iter, res, t0);
        if (res < tol) { record_iter_final(hist, iter, res, t0); return true; }
        x_old.swap(x);
    }
    record_iter_final(hist, iter > max_iter ? max_iter : iter, res, t0);
    return false;
}

static bool iterate_gs_sparse_core(const SparseCSR& A,
                                   const vector<double>& b,
                                   double tol,
                                   int max_iter,
                                   double omega,
                                   vector<double>& x,
                                   int& iter,
                                   double& res,
                                   IterHistory* hist)
{
    const int n = A.n_row;
    x.assign(n, 0.0);
    const auto t0 = chrono::high_resolution_clock::now();
    for (iter = 1; iter <= max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            const int diag_pos = (i < static_cast<int>(A.diag_pos.size()) ? A.diag_pos[i] : -1);
            if (diag_pos < 0 || fabs(A.val[diag_pos]) < 1.0e-30) return false;
            double s = b[i];
            for (int k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
                int j = A.col_id[k];
                if (j != i) s -= A.val[k] * x[j];
            }
            const double x_gs = s / A.val[diag_pos];
            x[i] = (1.0 - omega) * x[i] + omega * x_gs;
        }
        res = calculate_residual_sparse(A, x, b);
        record_iter_history(hist, iter, res, t0);
        if (res < tol) { record_iter_final(hist, iter, res, t0); return true; }
        if (!isfinite(res) || res > 1.0e80) return false;
    }
    record_iter_final(hist, iter > max_iter ? max_iter : iter, res, t0);
    return false;
}

static bool iterate_cg_sparse_core(const SparseCSR& A,
                                   const vector<double>& b,
                                   double tol,
                                   int max_iter,
                                   bool use_pcg,
                                   vector<double>& x,
                                   int& iter,
                                   double& res,
                                   IterHistory* hist)
{
    const int n = A.n_row;
    x.assign(n, 0.0);
    const auto t0 = chrono::high_resolution_clock::now();
    vector<double> r = b;
    vector<double> z(n, 0.0), p(n, 0.0), Ap(n, 0.0), diag(n, 1.0);

    if (use_pcg) {
        extract_diag_sparse(A, diag);
        for (int i = 0; i < n; ++i) {
            if (fabs(diag[i]) < 1.0e-30) diag[i] = 1.0;
            z[i] = r[i] / diag[i];
        }
        p = z;
    } else {
        p = r;
        z = r;
    }

    double rz_old = dot_vector(r, z);
    const double nb = max(1.0e-30, norm_vector(b));
    res = norm_vector(r) / nb;
    record_iter_history(hist, 0, res, t0);
    if (res < tol) { iter = 0; record_iter_final(hist, iter, res, t0); return true; }

    for (iter = 1; iter <= max_iter; ++iter) {
        A.multiply_vector(p, Ap);
        const double denom = dot_vector(p, Ap);
        if (fabs(denom) < 1.0e-30) return false;
        const double alpha = rz_old / denom;

        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        res = norm_vector(r) / nb;
        record_iter_history(hist, iter, res, t0);
        if (res < tol) { record_iter_final(hist, iter, res, t0); return true; }

        if (use_pcg) {
            for (int i = 0; i < n; ++i) z[i] = r[i] / diag[i];
        } else {
            z = r;
        }

        const double rz_new = dot_vector(r, z);
        const double beta = rz_new / max(1.0e-300, rz_old);
        for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        rz_old = rz_new;
    }
    record_iter_final(hist, iter > max_iter ? max_iter : iter, res, t0);
    return false;
}

SolverInfo solve_jacobi_dense(ModelData& model, const SolverParam& solver, TimeRecord& timer, const string& mode_name)
{
    const string step = "solve_jacobi_dense";
    timer.start_time(step);
    vector<double> A, b, x;
    int n = 0, iter = 0;
    double res = 1.0e100;
    prepare_dense_free(model, A, b, n);
    SolverInfo info = make_iter_info("jacobi_dense", "dense", 0.0, model.n_dof, n);
    IterHistory hist;
    init_iter_history(hist, info);
    info.converged = iterate_jacobi_dense_core(A, b, n, solver.tol, solver.max_iter, x, iter, res, &hist);
    info.iterations = iter;
    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    finish_iter(model, x, info, res);
    push_iter_history(model, hist);
    (void)mode_name;
    return info;
}

SolverInfo solve_gauss_seidel_dense(ModelData& model, const SolverParam& solver, TimeRecord& timer, const string& mode_name)
{
    const string step = "solve_gauss_seidel_dense";
    timer.start_time(step);
    vector<double> A, b, x;
    int n = 0, iter = 0;
    double res = 1.0e100;
    prepare_dense_free(model, A, b, n);
    SolverInfo info = make_iter_info("gauss_seidel_dense", "dense", 0.0, model.n_dof, n);
    IterHistory hist;
    init_iter_history(hist, info);
    info.converged = iterate_gs_dense_core(A, b, n, solver.tol, solver.max_iter, 1.0, x, iter, res, &hist);
    info.iterations = iter;
    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    finish_iter(model, x, info, res);
    push_iter_history(model, hist);
    (void)mode_name;
    return info;
}

SolverInfo solve_sor_dense(ModelData& model, const SolverParam& solver, double omega, TimeRecord& timer, const string& mode_name)
{
    const string step = "solve_sor_dense";
    timer.start_time(step);
    vector<double> A, b, x;
    int n = 0, iter = 0;
    double res = 1.0e100;
    prepare_dense_free(model, A, b, n);
    SolverInfo info = make_iter_info("sor_dense", "dense", omega, model.n_dof, n);
    IterHistory hist;
    init_iter_history(hist, info);
    info.converged = iterate_gs_dense_core(A, b, n, solver.tol, solver.max_iter, omega, x, iter, res, &hist);
    info.iterations = iter;
    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    finish_iter(model, x, info, res);
    push_iter_history(model, hist);
    (void)mode_name;
    return info;
}

SolverInfo solve_cg_dense(ModelData& model, const SolverParam& solver, TimeRecord& timer, const string& mode_name)
{
    const string step = "solve_cg_dense";
    timer.start_time(step);
    vector<double> A, b, x;
    int n = 0, iter = 0;
    double res = 1.0e100;
    prepare_dense_free(model, A, b, n);
    SolverInfo info = make_iter_info("cg_dense", "dense", 0.0, model.n_dof, n);
    IterHistory hist;
    init_iter_history(hist, info);
    info.converged = iterate_cg_dense_core(A, b, n, solver.tol, solver.max_iter, false, x, iter, res, &hist);
    info.iterations = iter;
    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    finish_iter(model, x, info, res);
    push_iter_history(model, hist);
    (void)mode_name;
    return info;
}

SolverInfo solve_pcg_dense(ModelData& model, const SolverParam& solver, TimeRecord& timer, const string& mode_name)
{
    const string step = "solve_pcg_dense";
    timer.start_time(step);
    vector<double> A, b, x;
    int n = 0, iter = 0;
    double res = 1.0e100;
    prepare_dense_free(model, A, b, n);
    SolverInfo info = make_iter_info("pcg_dense", "dense", 0.0, model.n_dof, n);
    IterHistory hist;
    init_iter_history(hist, info);
    info.converged = iterate_cg_dense_core(A, b, n, solver.tol, solver.max_iter, true, x, iter, res, &hist);
    info.iterations = iter;
    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    finish_iter(model, x, info, res);
    push_iter_history(model, hist);
    (void)mode_name;
    return info;
}

SolverInfo solve_jacobi_sparse(ModelData& model, const SolverParam& solver, TimeRecord& timer, const string& mode_name)
{
    const string step = "solve_jacobi_sparse";
    timer.start_time(step);
    vector<double> b, x;
    SparseCSR A = prepare_sparse_free(model, b);
    int iter = 0; double res = 1.0e100;
    SolverInfo info = make_iter_info("jacobi_sparse", "sparse", 0.0, model.n_dof, A.n_row);
    IterHistory hist;
    init_iter_history(hist, info);
    info.converged = iterate_jacobi_sparse_core(A, b, solver.tol, solver.max_iter, x, iter, res, &hist);
    info.iterations = iter;
    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    finish_iter(model, x, info, res);
    push_iter_history(model, hist);
    (void)mode_name;
    return info;
}

SolverInfo solve_gauss_seidel_sparse(ModelData& model, const SolverParam& solver, TimeRecord& timer, const string& mode_name)
{
    const string step = "solve_gauss_seidel_sparse";
    timer.start_time(step);
    vector<double> b, x;
    SparseCSR A = prepare_sparse_free(model, b);
    A.create_diag_pos();
    int iter = 0; double res = 1.0e100;
    SolverInfo info = make_iter_info("gauss_seidel_sparse", "sparse", 0.0, model.n_dof, A.n_row);
    IterHistory hist;
    init_iter_history(hist, info);
    info.converged = iterate_gs_sparse_core(A, b, solver.tol, solver.max_iter, 1.0, x, iter, res, &hist);
    info.iterations = iter;
    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    finish_iter(model, x, info, res);
    push_iter_history(model, hist);
    (void)mode_name;
    return info;
}

SolverInfo solve_sor_sparse(ModelData& model, const SolverParam& solver, double omega, TimeRecord& timer, const string& mode_name)
{
    const string step = "solve_sor_sparse";
    timer.start_time(step);
    vector<double> b, x;
    SparseCSR A = prepare_sparse_free(model, b);
    A.create_diag_pos();
    int iter = 0; double res = 1.0e100;
    SolverInfo info = make_iter_info("sor_sparse", "sparse", omega, model.n_dof, A.n_row);
    IterHistory hist;
    init_iter_history(hist, info);
    info.converged = iterate_gs_sparse_core(A, b, solver.tol, solver.max_iter, omega, x, iter, res, &hist);
    info.iterations = iter;
    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    finish_iter(model, x, info, res);
    push_iter_history(model, hist);
    (void)mode_name;
    return info;
}

SolverInfo solve_cg_sparse(ModelData& model, const SolverParam& solver, TimeRecord& timer, const string& mode_name)
{
    const string step = "solve_cg_sparse";
    timer.start_time(step);
    vector<double> b, x;
    SparseCSR A = prepare_sparse_free(model, b);
    int iter = 0; double res = 1.0e100;
    SolverInfo info = make_iter_info("cg_sparse", "sparse", 0.0, model.n_dof, A.n_row);
    IterHistory hist;
    init_iter_history(hist, info);
    info.converged = iterate_cg_sparse_core(A, b, solver.tol, solver.max_iter, false, x, iter, res, &hist);
    info.iterations = iter;
    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    finish_iter(model, x, info, res);
    push_iter_history(model, hist);
    (void)mode_name;
    return info;
}

SolverInfo solve_pcg_sparse(ModelData& model, const SolverParam& solver, TimeRecord& timer, const string& mode_name)
{
    const string step = "solve_pcg_sparse";
    timer.start_time(step);
    vector<double> b, x;
    SparseCSR A = prepare_sparse_free(model, b);
    int iter = 0; double res = 1.0e100;
    SolverInfo info = make_iter_info("pcg_sparse", "sparse", 0.0, model.n_dof, A.n_row);
    IterHistory hist;
    init_iter_history(hist, info);
    info.converged = iterate_cg_sparse_core(A, b, solver.tol, solver.max_iter, true, x, iter, res, &hist);
    info.iterations = iter;
    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    finish_iter(model, x, info, res);
    push_iter_history(model, hist);
    (void)mode_name;
    return info;
}

void solve_iterative(ModelData& model,
                     const SolverParam& solver,
                     TimeRecord& timer,
                     const string& mode_name)
{
    if (solver.run_jacobi_dense && !model.K_dense.empty()) solve_jacobi_dense(model, solver, timer, mode_name);
    if (solver.run_gauss_seidel_dense && !model.K_dense.empty()) solve_gauss_seidel_dense(model, solver, timer, mode_name);
    if (solver.run_sor_dense && !model.K_dense.empty()) {
        for (size_t i = 0; i < solver.omegas.size(); ++i) solve_sor_dense(model, solver, solver.omegas[i], timer, mode_name);
    }
    if (solver.run_cg_dense && !model.K_dense.empty()) solve_cg_dense(model, solver, timer, mode_name);
    if (solver.run_pcg_dense && !model.K_dense.empty()) solve_pcg_dense(model, solver, timer, mode_name);

    if (model.K_csr.n_row > 0) {
        if (solver.run_jacobi_sparse) solve_jacobi_sparse(model, solver, timer, mode_name);
        if (solver.run_gauss_seidel_sparse) solve_gauss_seidel_sparse(model, solver, timer, mode_name);
        if (solver.run_sor_sparse) {
            for (size_t i = 0; i < solver.omegas.size(); ++i) solve_sor_sparse(model, solver, solver.omegas[i], timer, mode_name);
        }
        if (solver.run_cg_sparse) solve_cg_sparse(model, solver, timer, mode_name);
        if (solver.run_pcg_sparse) solve_pcg_sparse(model, solver, timer, mode_name);
    }
}
