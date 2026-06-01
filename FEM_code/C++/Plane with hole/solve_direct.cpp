#include "fem_all.h"

#include <iostream>
#include <algorithm>

using namespace std;

static SolverInfo make_direct_info(const string& method,
                                   const string& matrix_type,
                                   int n_dof,
                                   int n_solve_dof)
{
    SolverInfo info;
    info.method = method;
    info.matrix_type = matrix_type;
    info.omega = 0.0;
    info.iterations = 1;
    info.n_dof = n_dof;
    info.n_solve_dof = n_solve_dof;
    return info;
}

SolverInfo solve_direct_dense_full(ModelData& model,
                                   TimeRecord& timer,
                                   const string& mode_name)
{
    const string step = "solve_direct_dense_full";
    timer.start_time(step);

    SolverInfo info = make_direct_info("direct_dense_full", "dense", model.n_dof, model.n_dof);

    vector<double> U;
    bool ok = solve_gauss_dense(model.K_dense, model.F, U, model.n_dof);
    info.converged = ok;

    if (ok) {
        model.U = U;
        model.U_ref = U;
        model.solve_success = true;
        info.final_residual = calculate_residual_dense(model.K_dense, U, model.F, model.n_dof);
        info.max_displacement = calculate_max_displacement(U);
    } else {
        info.final_residual = 1.0e100;
        model.add_warning("direct_dense_full failed");
    }

    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    model.direct_solver_info.push_back(info);
    model.all_solver_info.push_back(info);
    (void)mode_name;
    return info;
}

SolverInfo solve_direct_dense_block(ModelData& model,
                                    TimeRecord& timer,
                                    const string& mode_name)
{
    const string step = "solve_direct_dense_block";
    timer.start_time(step);

    const int nf = static_cast<int>(model.free_dof.size());
    SolverInfo info = make_direct_info("direct_dense_block", "dense", model.n_dof, nf);

    vector<double> Kff, Ff, Uf;
    extract_free_dof_dense(model.K_dense, model.F, model.free_dof, model.n_dof, Kff, Ff);

    bool ok = solve_gauss_dense(Kff, Ff, Uf, nf);
    info.converged = ok;

    if (ok) {
        recover_full_vector(Uf, model.free_dof, model.n_dof, model.U);
        model.U_ref = model.U;
        model.solve_success = true;
        info.final_residual = calculate_residual_dense(Kff, Uf, Ff, nf);
        info.max_displacement = calculate_max_displacement(model.U);
    } else {
        info.final_residual = 1.0e100;
        model.add_warning("direct_dense_block failed");
    }

    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    model.direct_solver_info.push_back(info);
    model.all_solver_info.push_back(info);
    (void)mode_name;
    return info;
}

SolverInfo solve_direct_sparse_block_to_dense(ModelData& model,
                                              TimeRecord& timer,
                                              const string& mode_name)
{
    const string step = "solve_direct_sparse_block_to_dense";
    timer.start_time(step);

    const int nf = static_cast<int>(model.free_dof.size());
    SolverInfo info = make_direct_info("direct_sparse_block_to_dense", "sparse", model.n_dof, nf);

    SparseCSR Kff;
    vector<double> Ff, Kff_dense, Uf;
    extract_free_dof_sparse(model.K_csr, model.F, model.free_dof, Kff, Ff);
    convert_csr_dense(Kff, Kff_dense);

    bool ok = solve_gauss_dense(Kff_dense, Ff, Uf, nf);
    info.converged = ok;

    if (ok) {
        recover_full_vector(Uf, model.free_dof, model.n_dof, model.U);
        if (model.U_ref.empty() || calculate_max_displacement(model.U_ref) == 0.0) model.U_ref = model.U;
        model.solve_success = true;
        info.final_residual = calculate_residual_dense(Kff_dense, Uf, Ff, nf);
        info.max_displacement = calculate_max_displacement(model.U);
    } else {
        info.final_residual = 1.0e100;
        model.add_warning("direct_sparse_block_to_dense failed");
    }

    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    model.direct_solver_info.push_back(info);
    model.all_solver_info.push_back(info);
    (void)mode_name;
    return info;
}

SolverInfo solve_cholesky_dense(ModelData& model,
                                TimeRecord& timer,
                                const string& mode_name)
{
    const string step = "solve_cholesky_dense";
    timer.start_time(step);

    const int nf = static_cast<int>(model.free_dof.size());
    SolverInfo info = make_direct_info("direct_cholesky_dense", "dense", model.n_dof, nf);

    vector<double> Kff, Ff, Uf;
    extract_free_dof_dense(model.K_dense, model.F, model.free_dof, model.n_dof, Kff, Ff);

    bool ok = solve_cholesky_dense_core(Kff, Ff, Uf, nf);
    info.converged = ok;

    if (ok) {
        recover_full_vector(Uf, model.free_dof, model.n_dof, model.U);
        if (model.U_ref.empty() || calculate_max_displacement(model.U_ref) == 0.0) model.U_ref = model.U;
        model.solve_success = true;
        info.final_residual = calculate_residual_dense(Kff, Uf, Ff, nf);
        info.max_displacement = calculate_max_displacement(model.U);
    } else {
        info.final_residual = 1.0e100;
        model.add_warning("Cholesky skipped or failed: matrix may not be positive definite");
    }

    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    model.direct_solver_info.push_back(info);
    model.all_solver_info.push_back(info);
    (void)mode_name;
    return info;
}

SolverInfo solve_ldlt_dense(ModelData& model,
                            TimeRecord& timer,
                            const string& mode_name)
{
    const string step = "solve_ldlt_dense";
    timer.start_time(step);

    SolverInfo info = make_direct_info("direct_ldlt_dense", "dense", model.n_dof,
                                       static_cast<int>(model.free_dof.size()));
    info.converged = false;
    info.iterations = 0;
    info.final_residual = 1.0e100;
    model.add_warning("LDLT interface is reserved but not enabled in this version");

    timer.stop_time(step);
    info.time_seconds = timer.get_time(step);
    model.direct_solver_info.push_back(info);
    model.all_solver_info.push_back(info);
    (void)mode_name;
    return info;
}

void solve_direct(ModelData& model,
                  const SolverParam& solver,
                  TimeRecord& timer,
                  const string& mode_name)
{
    if (solver.run_direct_dense_full && !model.K_dense.empty()) {
        solve_direct_dense_full(model, timer, mode_name);
    }
    if (solver.run_direct_dense_block && !model.K_dense.empty()) {
        solve_direct_dense_block(model, timer, mode_name);
    }
    if (solver.run_direct_sparse_block_to_dense && model.K_csr.n_row > 0) {
        solve_direct_sparse_block_to_dense(model, timer, mode_name);
    }
    if (solver.run_cholesky_dense && !model.K_dense.empty()) {
        solve_cholesky_dense(model, timer, mode_name);
    }
    if (solver.run_ldlt_dense && !model.K_dense.empty()) {
        solve_ldlt_dense(model, timer, mode_name);
    }
}
