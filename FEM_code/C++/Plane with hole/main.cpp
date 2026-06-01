#include "fem_all.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

// ============================================================
// main.cpp
// ------------------------------------------------------------
// 功能：
// 1. 统一管理全部可调参数。
// 2. 循环计算不同椭圆孔轴比 b/a。
// 3. 调用网格、单刚、总刚、边界、求解、后处理、输出模块。
// 4. 输出单工况结果和全部工况汇总结果。
// ============================================================

//生成输出文件路径
static string join_path_main(const string& dir, const string& name)
{
    if (dir.empty()) return name;
    char last = dir[dir.size() - 1];
    if (last == '/' || last == '\\') return dir + name;
    return dir + "/" + name;
}

//根据椭圆轴比ratio生成对应输出文件名称
static string make_ratio_folder(double ratio)
{
    ostringstream ss;
    ss << "ratio_" << fixed << setprecision(0) << ratio;
    return ss.str();
}

//搜索关键词, 输出不同算法运行的时间
static double sum_time_by_keyword(const TimeRecord& timer, const string& keyword)
{
    double total = 0.0;
    const map<string, double>& table = timer.get_table();

    for (const auto& item : table) {
        if (item.first.find(keyword) != string::npos) {
            total += item.second;
        }
    }

    return total;
}

// 生成对比结果
static CompareBaseFast create_compare_base_fast(const ModelData& model,
                                                const TimeRecord& timer,
                                                double ratio)
{
    CompareBaseFast info;
    info.ratio = ratio;
    info.n_node = model.n_node;
    info.n_elem = model.n_elem;

    info.error_U = 0.0;
    info.error_stress = 0.0;
    info.error_Kt_vm = 0.0;

    info.time_base = sum_time_by_keyword(timer, "base");
    info.time_fast = sum_time_by_keyword(timer, "fast");
    if (info.time_fast > 1.0e-14) info.speedup = info.time_base / info.time_fast;
    else info.speedup = 0.0;

    return info;
}

// 生成求解器结果对比信息
static void create_compare_solver_result(const ModelData& model,
                                         double ratio,
                                         vector<CompareSolverResult>& compare_infos)
{
    if (model.all_solver_info.empty()) return;

    int ref_id = -1;
    for (size_t i = 0; i < model.all_solver_info.size(); ++i) {
        if (model.all_solver_info[i].converged) {
            ref_id = static_cast<int>(i);
            break;
        }
    }
    if (ref_id < 0) return;

    const SolverInfo& ref = model.all_solver_info[static_cast<size_t>(ref_id)];
    const double ref_disp = max(1.0e-30, fabs(ref.max_displacement));
    const double ref_res = max(1.0e-30, fabs(ref.final_residual));

    for (size_t i = 0; i < model.all_solver_info.size(); ++i) {
        const SolverInfo& cur = model.all_solver_info[i];
        CompareSolverResult c;
        c.ratio = ratio;
        c.method_ref = ref.method;
        c.method_compare = cur.method;

        // 由于 SolverInfo 中不保存每个求解器的完整位移向量，
        // 这里的 error_U 先写 0，主要对比最大位移和残差。
        c.error_U = 0.0;
        c.error_max_disp = fabs(cur.max_displacement - ref.max_displacement) / ref_disp;
        c.error_residual = fabs(cur.final_residual - ref.final_residual) / ref_res;
        compare_infos.push_back(c);
    }
}
// 统计总耗时
static void append_time_summary(const TimeRecord& timer,
                                double ratio,
                                vector<pair<string, double> >& time_infos)
{
    const map<string, double>& table = timer.get_table();
    for (map<string, double>::const_iterator it = table.begin(); it != table.end(); ++it) {
        ostringstream name;
        name << "ratio_" << fixed << setprecision(0) << ratio << "_" << it->first;
        time_infos.push_back(make_pair(name.str(), it->second));
    }
}

// 防止卡死
static void apply_large_problem_safety(SolverParam& solver,
                                       int n_dof_est,
                                       bool auto_skip_expensive_solver,
                                       int expensive_dof_limit)
{
    if (!auto_skip_expensive_solver) return;
    if (n_dof_est <= expensive_dof_limit) return;

    cout << "[Safety] n_dof = " << n_dof_est
         << " is larger than expensive_dof_limit = " << expensive_dof_limit << endl;
    cout << "[Safety] Skip dense_full/Jacobi/GS/SOR for practical running. "
         << "Set auto_skip_expensive_solver=false to force all solvers." << endl;

    // 完整 dense 高斯消元、dense 子矩阵直接法和低效迭代法对大规模有限元矩阵非常慢。
    // 默认保留 sparse CG / sparse PCG 作为大网格主要求解器。
    solver.run_direct_dense_full = false;
    solver.run_direct_dense_block = false;
    solver.run_direct_sparse_block_to_dense = false;
    solver.run_cholesky_dense = false;
    solver.run_ldlt_dense = false;

    solver.run_jacobi_dense = false;
    solver.run_gauss_seidel_dense = false;
    solver.run_sor_dense = false;
    solver.run_cg_dense = false;
    solver.run_pcg_dense = false;

    solver.run_jacobi_sparse = false;
    solver.run_gauss_seidel_sparse = false;
    solver.run_sor_sparse = false;
}

int main()
{
    cout << "============================================================" << endl;
    cout << "FEM Q4 Half Plate With Center Ellipse Hole" << endl;
    cout << "二维平面应力：带中心椭圆孔矩形板 1/2 模型" << endl;
    cout << "============================================================" << endl;

    // ========================================================
    // 一、可调参数
    // ========================================================

    // ---------------- 几何参数 ----------------
    double L = 84.0;       // 板长，单位 mm
    double H = 84.0;       // 完整板高度，单位 mm；本程序计算上半板 H/2
    double R0 = 6.0;       // 面积等效圆孔半径，单位 mm

    // ---------------- 材料参数 ----------------
    double E = 110000.0;   // 杨氏模量，单位 MPa = N/mm^2
    double nu = 0.3;       // 泊松比
    double t = 1.0;        // 厚度，单位 mm

    // ---------------- 载荷参数 ----------------
    double q = 1.0;        // 右端均布拉伸载荷，单位 MPa = N/mm^2

    // ---------------- 网格参数 ----------------
    int n_theta = 64;      // TFI 周向总分段数，必须能被 4 整除
    int n_inner = 8;       // TFI 径向分段的一部分；与 n_outer 相加得到 n_eta
    int n_outer = 16;      // TFI 径向分段的一部分；n_eta = n_inner + n_outer，推荐总数约 24
    double lambda = 1.8;   // TFI 网格中保留该参数用于兼容旧接口，实际不再控制外相似椭圆

    // ---------------- 求解参数 ----------------
    double tol = 1.0e-8;
    int max_iter = 200000;

    // ---------------- 椭圆孔轴比工况 ----------------
    vector<double> ratios;
    ratios.push_back(1.0);
    ratios.push_back(2.0);
    ratios.push_back(3.0);
    ratios.push_back(4.0);
    ratios.push_back(5.0);

    // ---------------- SOR 松弛因子 ----------------
    vector<double> omegas;
    for (int i = 10; i <= 19; ++i) omegas.push_back(static_cast<double>(i) / 10.0);

    // ---------------- 运行模式 ----------------
    RunMode run_mode = BOTH_MODE;
    MatrixMode matrix_mode = BOTH_MATRIX;
    AssembleMode assemble_mode = ASSEMBLE_ALL;
    StressMode stress_mode = CENTER_STRESS;

    // ---------------- 加速开关 ----------------
    bool use_sparse_matrix = true;
    bool use_fast_mesh = true;
    bool use_fast_assemble = true;
    bool use_fast_boundary = true;
    bool use_fast_post = true;
    bool use_fast_io = true;
    bool use_precompute_gauss = true;
    bool use_precompute_edof = true;
    bool use_static_local_array = true;
    bool use_parallel_loop = false;
    bool compare_base_fast = true;
    bool compare_solver_result = true;

    // ---------------- 求解器开关 ----------------
    bool run_direct_dense_full = true;
    bool run_direct_dense_block = true;
    bool run_direct_sparse_block_to_dense = true;
    bool run_cholesky_dense = true;
    bool run_ldlt_dense = false;

    bool run_jacobi_dense = true;
    bool run_gauss_seidel_dense = true;
    bool run_sor_dense = true;
    bool run_cg_dense = true;
    bool run_pcg_dense = true;

    bool run_jacobi_sparse = true;
    bool run_gauss_seidel_sparse = true;
    bool run_sor_sparse = true;
    bool run_cg_sparse = true;
    bool run_pcg_sparse = true;

 
    bool auto_skip_expensive_solver = false;   // TFI 正式网格自由度较大，默认跳过高耗时求解器，避免程序卡死
    int expensive_dof_limit = 1500;

    // ---------------- 计时精度 ----------------
    int timer_detail_level = 2;

    // ---------------- 输出路径 ----------------
    string output_root = "output";

    create_folder(output_root);

    // ========================================================
    // 二、全部工况汇总容器
    // ========================================================

    vector<CaseSummary> all_case_summary;
    vector<SolverInfo> all_solver_summary;
    vector<pair<string, double> > all_time_summary;
    vector<CompareBaseFast> all_compare_base_fast;
    vector<CompareSolverResult> all_compare_solver_result;

    // ========================================================
    // 三、循环计算每个 b/a 工况
    // ========================================================

    for (size_t ir = 0; ir < ratios.size(); ++ir) {
        double ratio = ratios[ir];
        cout << "\n============================================================" << endl;
        cout << "Start ratio = " << ratio << endl;
        cout << "============================================================" << endl;

        TimeRecord timer;
        timer.start_time("total_case_time");

        // ---------------- 1. 创建材料参数 ----------------
        timer.start_time("build_material");
        MaterialParam mat(E, nu, t);
        timer.stop_time("build_material");

        // ---------------- 2. 创建几何参数 ----------------
        timer.start_time("build_geometry");
        GeometryParam geom(L, H, R0, ratio);
        timer.stop_time("build_geometry");

        // ---------------- 3. 创建网格参数 ----------------
        MeshParam mesh(n_theta, n_inner, n_outer, lambda, use_fast_mesh);

        // ---------------- 4. 创建载荷参数 ----------------
        LoadParam load(q, 0);

        // ---------------- 5. 创建求解器参数 ----------------
        SolverParam solver;
        solver.tol = tol;
        solver.max_iter = max_iter;
        solver.omegas = omegas;

        solver.run_direct_dense_full = run_direct_dense_full;
        solver.run_direct_dense_block = run_direct_dense_block;
        solver.run_direct_sparse_block_to_dense = run_direct_sparse_block_to_dense;
        solver.run_cholesky_dense = run_cholesky_dense;
        solver.run_ldlt_dense = run_ldlt_dense;

        solver.run_jacobi_dense = run_jacobi_dense;
        solver.run_gauss_seidel_dense = run_gauss_seidel_dense;
        solver.run_sor_dense = run_sor_dense;
        solver.run_cg_dense = run_cg_dense;
        solver.run_pcg_dense = run_pcg_dense;

        solver.run_jacobi_sparse = run_jacobi_sparse;
        solver.run_gauss_seidel_sparse = run_gauss_seidel_sparse;
        solver.run_sor_sparse = run_sor_sparse;
        solver.run_cg_sparse = run_cg_sparse;
        solver.run_pcg_sparse = run_pcg_sparse;

        // ---------------- 6. 创建加速参数 ----------------
        SpeedParam speed;
        speed.run_mode = run_mode;
        speed.matrix_mode = matrix_mode;
        speed.assemble_mode = assemble_mode;
        speed.stress_mode = stress_mode;
        speed.use_sparse_matrix = use_sparse_matrix;
        speed.use_fast_mesh = use_fast_mesh;
        speed.use_fast_assemble = use_fast_assemble;
        speed.use_fast_boundary = use_fast_boundary;
        speed.use_fast_post = use_fast_post;
        speed.use_fast_io = use_fast_io;
        speed.use_precompute_gauss = use_precompute_gauss;
        speed.use_precompute_edof = use_precompute_edof;
        speed.use_static_local_array = use_static_local_array;
        speed.use_parallel_loop = use_parallel_loop;
        speed.compare_base_fast = compare_base_fast;
        speed.compare_solver_result = compare_solver_result;
        speed.timer_detail_level = timer_detail_level;

        // 打印参数，方便截图提交
        mat.print_material();
        geom.print_geometry();
        mesh.print_mesh_param();
        load.print_load();
        solver.print_solver_param();
        speed.print_speed_param();

        // ---------------- 7. 生成网格 ----------------
        ModelData model;
        create_mesh(model, geom, mesh, speed, timer, "total");

        // 根据网格规模自动保护，避免极慢求解器拖死调试过程
        apply_large_problem_safety(solver, model.n_dof,
                                   auto_skip_expensive_solver,
                                   expensive_dof_limit);

        // ---------------- 8. 检查网格质量 ----------------
        check_mesh(model, mesh, timer, "total");
        model.print_model_size();
        print_mesh_info(model);

        // ---------------- 9. 预计算高斯点 ----------------
        create_Gauss(model.gauss_data, timer, "total");

        // ---------------- 10. 预计算单元自由度编号 ----------------
        if (speed.use_precompute_edof) {
            create_edof(model, timer, "total");
        }

        // ---------------- 11. 组装总体刚度矩阵 ----------------
        create_K_global(model, mat, speed, timer, "total");

        // ---------------- 12. 施加右端均布载荷 ----------------
        apply_load(model, geom, load, timer, "total");
        check_total_load(model, geom, mat, load);

        // ---------------- 13-14. 建立自由度集合并施加边界条件 ----------------
        apply_boundary(model, speed, timer, "total");

        // ---------------- 15. 直接法求解 ----------------
        solve_direct(model, solver, timer, "total");

        // ---------------- 16. 迭代法求解 ----------------
        solve_iterative(model, solver, timer, "total");

        // ---------------- 17. 选择参考解 ----------------
        if (model.U_ref.empty() || calculate_max_displacement(model.U_ref) <= 0.0) {
            model.U_ref = model.U;
        }

        // ---------------- 18. 后处理计算应变、应力 ----------------
        calculate_result(model, mat, load, speed, timer, "total");

        // ---------------- 19. 工况汇总 ----------------
        calculate_summary(model, geom, load);

        // ---------------- 20. 求解器结果对比 ----------------
        if (speed.compare_solver_result) {
            create_compare_solver_result(model, ratio, all_compare_solver_result);
        }

        // ---------------- 21. 输出单工况结果 ----------------
        string case_dir = join_path_main(output_root, make_ratio_folder(ratio));
        create_folder(case_dir);
        write_result(model, geom, mat, load, speed, timer, case_dir);

        // ---------------- 22. 结束单工况计时并再次写入 time_summary.csv ----------------
        timer.stop_time("total_case_time");
        timer.write_time(join_path_main(case_dir, "time_summary.csv"), "total");
        timer.print_time();

        // ---------------- 23. 保存全局汇总数据 ----------------
        all_case_summary.push_back(model.summary);
        for (size_t i = 0; i < model.all_solver_info.size(); ++i) {
            SolverInfo tmp = model.all_solver_info[i];
            all_solver_summary.push_back(tmp);
        }
        append_time_summary(timer, ratio, all_time_summary);

        if (speed.compare_base_fast) {
            all_compare_base_fast.push_back(create_compare_base_fast(model, timer, ratio));
        }

        cout << "Finish ratio = " << ratio << endl;
    }

    // ========================================================
    // 四、输出全部工况汇总文件
    // ========================================================

    write_summary_ratio(all_case_summary, join_path_main(output_root, "summary_ratio.csv"));
    write_summary_solver(all_solver_summary, join_path_main(output_root, "summary_solver.csv"));
    write_summary_time(all_time_summary, join_path_main(output_root, "summary_time.csv"));
    write_compare_base_fast(all_compare_base_fast, join_path_main(output_root, "compare_base_fast.csv"));
    write_compare_solver_result(all_compare_solver_result, join_path_main(output_root, "compare_solver_result.csv"));

    cout << "\n============================================================" << endl;
    cout << "All cases finished." << endl;
    cout << "Results are written to folder: " << output_root << endl;
    cout << "============================================================" << endl;

    return 0;
}
