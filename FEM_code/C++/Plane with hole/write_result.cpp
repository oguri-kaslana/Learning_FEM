#include "fem_all.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cerrno>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#endif

using namespace std;

bool create_folder(const string& folder)
{
    if (folder.empty()) return false;

#ifdef _WIN32
    int ret = _mkdir(folder.c_str());
#else
    int ret = mkdir(folder.c_str(), 0755);
#endif
    if (ret == 0 || errno == EEXIST) return true;

    // 支持简单的多级目录：逐段创建。
    string current;
    for (size_t i = 0; i < folder.size(); ++i) {
        char c = folder[i];
        current.push_back(c);
        if (c == '/' || c == '\\' || i + 1 == folder.size()) {
            if (current.size() <= 1) continue;
#ifdef _WIN32
            _mkdir(current.c_str());
#else
            mkdir(current.c_str(), 0755);
#endif
        }
    }
    return true;
}

static string join_path(const string& dir, const string& name)
{
    if (dir.empty()) return name;
    const char last = dir[dir.size() - 1];
    if (last == '/' || last == '\\') return dir + name;
    return dir + "/" + name;
}

static void open_check(ofstream& fout, const string& filename)
{
    if (!fout.is_open()) {
        cerr << "[write_result warning] cannot open file: " << filename << endl;
    }
}

void write_nodes(const ModelData& model, const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "node_id,x,y,flag_left,flag_right,flag_top,flag_sym,flag_hole\n";
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        const NodeInfo& n = model.nodes[i];
        fout << n.id << ',' << n.x << ',' << n.y << ','
             << n.flag_left << ',' << n.flag_right << ',' << n.flag_top << ','
             << n.flag_sym << ',' << n.flag_hole << '\n';
    }
}

void write_elements(const ModelData& model, const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << "elem_id,n1,n2,n3,n4\n";
    for (size_t i = 0; i < model.elements.size(); ++i) {
        const ElementInfo& e = model.elements[i];
        fout << e.id << ',' << e.node[0] << ',' << e.node[1] << ','
             << e.node[2] << ',' << e.node[3] << '\n';
    }
}

void write_displacement(const ModelData& model, const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "node_id,x,y,ux,uy,umag\n";
    const vector<double>& U = model.U_ref.empty() ? model.U : model.U_ref;
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        const NodeInfo& n = model.nodes[i];
        const double ux = (2 * n.id < static_cast<int>(U.size())) ? U[2 * n.id] : 0.0;
        const double uy = (2 * n.id + 1 < static_cast<int>(U.size())) ? U[2 * n.id + 1] : 0.0;
        fout << n.id << ',' << n.x << ',' << n.y << ','
             << ux << ',' << uy << ',' << calculate_umag(ux, uy) << '\n';
    }
}

void write_strain(const ModelData& model, const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "elem_id,x_center,y_center,epsilon_x,epsilon_y,gamma_xy\n";
    for (size_t i = 0; i < model.elem_result.size(); ++i) {
        const ElementResult& r = model.elem_result[i];
        fout << r.elem_id << ',' << r.x_center << ',' << r.y_center << ','
             << r.epsilon_x << ',' << r.epsilon_y << ',' << r.gamma_xy << '\n';
    }
}

void write_stress(const ModelData& model, const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "elem_id,x_center,y_center,sigma_x,sigma_y,tau_xy,sigma_vm\n";
    for (size_t i = 0; i < model.elem_result.size(); ++i) {
        const ElementResult& r = model.elem_result[i];
        fout << r.elem_id << ',' << r.x_center << ',' << r.y_center << ','
             << r.sigma_x << ',' << r.sigma_y << ',' << r.tau_xy << ',' << r.sigma_vm << '\n';
    }
}

void write_case_summary(const ModelData& model, const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "ratio,a,b,max_ux,max_uy,max_umag,max_sigma_x,max_sigma_y,max_tau_xy,max_sigma_vm,Kt_x,Kt_vm\n";
    const CaseSummary& s = model.summary;
    fout << s.ratio << ',' << s.a << ',' << s.b << ','
         << s.max_ux << ',' << s.max_uy << ',' << s.max_umag << ','
         << s.max_sigma_x << ',' << s.max_sigma_y << ',' << s.max_tau_xy << ','
         << s.max_sigma_vm << ',' << s.Kt_x << ',' << s.Kt_vm << '\n';
}

void write_direct_solver_summary(const ModelData& model, const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "method,matrix_type,n_dof,n_solve_dof,time_seconds,final_residual,max_displacement\n";
    for (size_t i = 0; i < model.direct_solver_info.size(); ++i) {
        const SolverInfo& s = model.direct_solver_info[i];
        fout << s.method << ',' << s.matrix_type << ',' << s.n_dof << ',' << s.n_solve_dof << ','
             << s.time_seconds << ',' << s.final_residual << ',' << s.max_displacement << '\n';
    }
}

void write_iterative_solver_summary(const ModelData& model, const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "method,matrix_type,omega,converged,iterations,time_seconds,final_residual,max_displacement\n";
    for (size_t i = 0; i < model.iterative_solver_info.size(); ++i) {
        const SolverInfo& s = model.iterative_solver_info[i];
        fout << s.method << ',' << s.matrix_type << ',' << s.omega << ','
             << (s.converged ? 1 : 0) << ',' << s.iterations << ','
             << s.time_seconds << ',' << s.final_residual << ',' << s.max_displacement << '\n';
    }
}

void write_iterative_history(const ModelData& model, const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;

    fout << setprecision(12) << fixed;
    fout << "method,matrix_type,omega,iteration,residual,time_seconds\n";

    for (size_t h = 0; h < model.iterative_history.size(); ++h) {
        const IterHistory& hist = model.iterative_history[h];
        const size_t n = hist.iteration.size();
        for (size_t i = 0; i < n; ++i) {
            fout << hist.method << ','
                 << hist.matrix_type << ','
                 << hist.omega << ','
                 << hist.iteration[i] << ','
                 << hist.residual[i] << ','
                 << hist.time_seconds[i] << '\n';
        }
    }
}

void write_check_info(const ModelData& model,
                      const GeometryParam& geom,
                      const MaterialParam& mat,
                      const LoadParam& load,
                      const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "L = " << geom.L << "\n";
    fout << "H = " << geom.H << "\n";
    fout << "a = " << geom.a << "\n";
    fout << "b = " << geom.b << "\n";
    fout << "ratio = " << geom.ratio << "\n";
    fout << "E = " << mat.E << "\n";
    fout << "nu = " << mat.nu << "\n";
    fout << "t = " << mat.t << "\n";
    fout << "q = " << load.q << "\n";
    fout << "n_node = " << model.n_node << "\n";
    fout << "n_elem = " << model.n_elem << "\n";
    fout << "n_dof = " << model.n_dof << "\n";
    fout << "left_fixed_dof = " << model.n_left_fixed_dof << "\n";
    fout << "symmetry_fixed_dof = " << model.n_sym_fixed_dof << "\n";
    fout << "right_load_edge = " << model.n_right_load_edge << "\n";
    fout << "actual_total_force = " << model.total_load_actual << "\n";
    fout << "theory_total_force = " << model.total_load_theory << "\n";
    fout << "min_detJ = " << model.min_detJ << "\n";
    fout << "max_detJ = " << model.max_detJ << "\n";
    fout << "has_negative_detJ = " << (model.has_negative_detJ ? 1 : 0) << "\n";
    fout << "dense_matrix_size = " << model.n_dof << " x " << model.n_dof << "\n";
    fout << "sparse_nnz = " << model.sparse_nnz << "\n";
    fout << "sparse_rate = " << model.sparse_rate << "\n";
    fout << "solve_success = " << (model.solve_success ? 1 : 0) << "\n";
    fout << "warning_count = " << model.warning_info.size() << "\n";
    for (size_t i = 0; i < model.warning_info.size(); ++i) {
        fout << "warning_" << i + 1 << " = " << model.warning_info[i] << "\n";
    }
}

void write_result_base(const ModelData& model,
                       const GeometryParam& geom,
                       const MaterialParam& mat,
                       const LoadParam& load,
                       TimeRecord& timer,
                       const string& output_dir)
{
    timer.start_time("write_result_base");
    create_folder(output_dir);
    write_nodes(model, join_path(output_dir, "nodes.csv"));
    write_elements(model, join_path(output_dir, "elements.csv"));
    write_displacement(model, join_path(output_dir, "displacement.csv"));
    write_strain(model, join_path(output_dir, "strain.csv"));
    write_stress(model, join_path(output_dir, "stress.csv"));
    write_direct_solver_summary(model, join_path(output_dir, "direct_solver_summary.csv"));
    write_iterative_solver_summary(model, join_path(output_dir, "iterative_solver_summary.csv"));
    write_iterative_history(model, join_path(output_dir, "iterative_history.csv"));
    write_case_summary(model, join_path(output_dir, "case_summary.csv"));
    write_check_info(model, geom, mat, load, join_path(output_dir, "check_info.txt"));
    timer.write_time(join_path(output_dir, "time_summary.csv"), "base");
    timer.stop_time("write_result_base");
}

void write_result_fast(const ModelData& model,
                       const GeometryParam& geom,
                       const MaterialParam& mat,
                       const LoadParam& load,
                       TimeRecord& timer,
                       const string& output_dir)
{
    timer.start_time("write_result_fast");
    create_folder(output_dir);

    // 文件规模较大时仍调用分文件写入，内部统一使用 '\n' 和固定精度，减少 endl 刷新开销。
    write_nodes(model, join_path(output_dir, "nodes.csv"));
    write_elements(model, join_path(output_dir, "elements.csv"));
    write_displacement(model, join_path(output_dir, "displacement.csv"));
    write_strain(model, join_path(output_dir, "strain.csv"));
    write_stress(model, join_path(output_dir, "stress.csv"));
    write_direct_solver_summary(model, join_path(output_dir, "direct_solver_summary.csv"));
    write_iterative_solver_summary(model, join_path(output_dir, "iterative_solver_summary.csv"));
    write_iterative_history(model, join_path(output_dir, "iterative_history.csv"));
    write_case_summary(model, join_path(output_dir, "case_summary.csv"));
    write_check_info(model, geom, mat, load, join_path(output_dir, "check_info.txt"));
    timer.write_time(join_path(output_dir, "time_summary.csv"), "fast");

    timer.stop_time("write_result_fast");
}

void write_result(const ModelData& model,
                  const GeometryParam& geom,
                  const MaterialParam& mat,
                  const LoadParam& load,
                  const SpeedParam& speed,
                  TimeRecord& timer,
                  const string& output_dir)
{
    timer.start_time("write_result_total");
    if (speed.use_fast_io) write_result_fast(model, geom, mat, load, timer, output_dir);
    else write_result_base(model, geom, mat, load, timer, output_dir);
    timer.stop_time("write_result_total");
}

void write_summary_ratio(const vector<CaseSummary>& summaries,
                         const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "ratio,a,b,max_ux,max_uy,max_umag,max_sigma_x,max_sigma_y,max_tau_xy,max_sigma_vm,max_sigma_vm_elem,max_sigma_vm_x,max_sigma_vm_y,Kt_x,Kt_vm\n";
    for (size_t i = 0; i < summaries.size(); ++i) {
        const CaseSummary& s = summaries[i];
        fout << s.ratio << ',' << s.a << ',' << s.b << ','
             << s.max_ux << ',' << s.max_uy << ',' << s.max_umag << ','
             << s.max_sigma_x << ',' << s.max_sigma_y << ',' << s.max_tau_xy << ','
             << s.max_sigma_vm << ',' << s.max_sigma_vm_elem << ','
             << s.max_sigma_vm_x << ',' << s.max_sigma_vm_y << ','
             << s.Kt_x << ',' << s.Kt_vm << '\n';
    }
}

void write_summary_solver(const vector<SolverInfo>& solver_infos,
                          const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "method,matrix_type,omega,converged,iterations,time_seconds,final_residual,max_displacement,n_dof,n_solve_dof\n";
    for (size_t i = 0; i < solver_infos.size(); ++i) {
        const SolverInfo& s = solver_infos[i];
        fout << s.method << ',' << s.matrix_type << ',' << s.omega << ','
             << (s.converged ? 1 : 0) << ',' << s.iterations << ','
             << s.time_seconds << ',' << s.final_residual << ',' << s.max_displacement << ','
             << s.n_dof << ',' << s.n_solve_dof << '\n';
    }
}

void write_summary_time(const vector<pair<string, double> >& time_infos,
                        const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "step_name,time_seconds\n";
    for (size_t i = 0; i < time_infos.size(); ++i) {
        fout << time_infos[i].first << ',' << time_infos[i].second << '\n';
    }
}

void write_compare_base_fast(const vector<CompareBaseFast>& compare_infos,
                             const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "ratio,n_node,n_elem,error_U,error_stress,error_Kt_vm,time_base,time_fast,speedup\n";
    for (size_t i = 0; i < compare_infos.size(); ++i) {
        const CompareBaseFast& c = compare_infos[i];
        fout << c.ratio << ',' << c.n_node << ',' << c.n_elem << ','
             << c.error_U << ',' << c.error_stress << ',' << c.error_Kt_vm << ','
             << c.time_base << ',' << c.time_fast << ',' << c.speedup << '\n';
    }
}

void write_compare_solver_result(const vector<CompareSolverResult>& compare_infos,
                                 const string& filename)
{
    ofstream fout(filename.c_str());
    open_check(fout, filename);
    if (!fout.is_open()) return;
    fout << setprecision(12) << fixed;
    fout << "ratio,method_ref,method_compare,error_U,error_max_disp,error_residual\n";
    for (size_t i = 0; i < compare_infos.size(); ++i) {
        const CompareSolverResult& c = compare_infos[i];
        fout << c.ratio << ',' << c.method_ref << ',' << c.method_compare << ','
             << c.error_U << ',' << c.error_max_disp << ',' << c.error_residual << '\n';
    }
}
