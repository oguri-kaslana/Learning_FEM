#include "fem_all.h"

#include <iomanip>
#include <algorithm>

using std::cout;
using std::endl;
using std::string;
using std::vector;

// ============================================================
// 1. MaterialParam
// ============================================================

MaterialParam::MaterialParam()
{
    E = 110000.0;
    nu = 0.3;
    t = 1.0;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            D[i][j] = 0.0;
        }
    }

    build_material();
}

MaterialParam::MaterialParam(double E_in, double nu_in, double t_in)
{
    E = E_in;
    nu = nu_in;
    t = t_in;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            D[i][j] = 0.0;
        }
    }

    build_material();
}

void MaterialParam::build_material()
{
    // 平面应力本构矩阵：
    // D = E / (1 - nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2]
    const double coef = E / (1.0 - nu * nu);

    D[0][0] = coef;
    D[0][1] = coef * nu;
    D[0][2] = 0.0;

    D[1][0] = coef * nu;
    D[1][1] = coef;
    D[1][2] = 0.0;

    D[2][0] = 0.0;
    D[2][1] = 0.0;
    D[2][2] = coef * (1.0 - nu) / 2.0;
}

void MaterialParam::print_material() const
{
    cout << "================ MaterialParam ================" << endl;
    cout << "E  = " << E << " MPa" << endl;
    cout << "nu = " << nu << endl;
    cout << "t  = " << t << " mm" << endl;
    cout << "D matrix:" << endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            cout << std::setw(18) << D[i][j] << " ";
        }
        cout << endl;
    }
}

// ============================================================
// 2. GeometryParam
// ============================================================

GeometryParam::GeometryParam()
{
    L = 84.0;
    H = 84.0;
    R0 = 6.0;
    ratio = 1.0;

    cx = 0.0;
    cy = 0.0;
    a = 0.0;
    b = 0.0;

    build_geometry();
}

GeometryParam::GeometryParam(double L_in, double H_in, double R0_in, double ratio_in)
{
    L = L_in;
    H = H_in;
    R0 = R0_in;
    ratio = ratio_in;

    cx = 0.0;
    cy = 0.0;
    a = 0.0;
    b = 0.0;

    build_geometry();
}

void GeometryParam::build_geometry()
{
    // 完整板中心孔位于 x = L/2, y = 0。
    // 上半板模型只保留 y >= 0 的区域。
    cx = L / 2.0;
    cy = 0.0;

    if (ratio <= 0.0) {
        cout << "Warning: ratio <= 0, reset ratio = 1.0" << endl;
        ratio = 1.0;
    }

    // 保持孔洞面积不变：pi*a*b = pi*R0^2，且 ratio=b/a。
    a = R0 / std::sqrt(ratio);
    b = R0 * std::sqrt(ratio);
}

void GeometryParam::print_geometry() const
{
    const double pi = std::acos(-1.0);
    cout << "================ GeometryParam ================" << endl;
    cout << "L      = " << L << " mm" << endl;
    cout << "H      = " << H << " mm" << endl;
    cout << "cx     = " << cx << " mm" << endl;
    cout << "cy     = " << cy << " mm" << endl;
    cout << "R0     = " << R0 << " mm" << endl;
    cout << "ratio  = " << ratio << endl;
    cout << "a      = " << a << " mm" << endl;
    cout << "b      = " << b << " mm" << endl;
    cout << "area   = " << pi * a * b << " mm^2" << endl;
}

// ============================================================
// 3. MeshParam
// ============================================================

MeshParam::MeshParam()
{
    n_theta = 96;
    n_inner = 8;
    n_outer = 12;
    lambda = 1.8;

    min_detJ = std::numeric_limits<double>::max();
    max_detJ = -std::numeric_limits<double>::max();

    use_fast_mesh = true;
}

MeshParam::MeshParam(int n_theta_in,
                     int n_inner_in,
                     int n_outer_in,
                     double lambda_in,
                     bool use_fast_mesh_in)
{
    n_theta = n_theta_in;
    n_inner = n_inner_in;
    n_outer = n_outer_in;
    lambda = lambda_in;
    use_fast_mesh = use_fast_mesh_in;

    if (n_theta < 2) {
        cout << "Warning: n_theta < 2, reset n_theta = 2" << endl;
        n_theta = 2;
    }
    if (n_theta % 2 != 0) {
        cout << "Warning: n_theta is not even, reset n_theta = n_theta + 1" << endl;
        n_theta += 1;
    }
    if (n_inner < 1) {
        cout << "Warning: n_inner < 1, reset n_inner = 1" << endl;
        n_inner = 1;
    }
    if (n_outer < 1) {
        cout << "Warning: n_outer < 1, reset n_outer = 1" << endl;
        n_outer = 1;
    }
    if (lambda <= 1.0) {
        cout << "Warning: lambda <= 1.0, reset lambda = 1.8" << endl;
        lambda = 1.8;
    }

    min_detJ = std::numeric_limits<double>::max();
    max_detJ = -std::numeric_limits<double>::max();
}

void MeshParam::print_mesh_param() const
{
    cout << "================ MeshParam ================" << endl;
    cout << "n_theta      = " << n_theta << endl;
    cout << "n_inner      = " << n_inner << endl;
    cout << "n_outer      = " << n_outer << endl;
    cout << "lambda       = " << lambda << endl;
    cout << "min_detJ     = " << min_detJ << endl;
    cout << "max_detJ     = " << max_detJ << endl;
    cout << "use_fast_mesh= " << (use_fast_mesh ? "true" : "false") << endl;
}

// ============================================================
// 4. LoadParam
// ============================================================

LoadParam::LoadParam()
{
    q = 1.0;
    direction = 0;
    total_force = 0.0;
}

LoadParam::LoadParam(double q_in, int direction_in)
{
    q = q_in;
    direction = direction_in;
    total_force = 0.0;

    if (direction != 0 && direction != 1) {
        cout << "Warning: invalid load direction, reset direction = 0" << endl;
        direction = 0;
    }
}

void LoadParam::print_load() const
{
    cout << "================ LoadParam ================" << endl;
    cout << "q           = " << q << " N/mm^2" << endl;
    cout << "direction   = " << (direction == 0 ? "+x" : "+y") << endl;
    cout << "total_force = " << total_force << " N" << endl;
}

// ============================================================
// 5. SolverParam
// ============================================================

SolverParam::SolverParam()
{
    tol = 1.0e-8;
    max_iter = 200000;

    run_direct_dense_full = true;
    run_direct_dense_block = true;
    run_direct_sparse_block_to_dense = true;
    run_cholesky_dense = true;
    run_ldlt_dense = false;

    run_jacobi_dense = true;
    run_gauss_seidel_dense = true;
    run_sor_dense = true;
    run_cg_dense = true;
    run_pcg_dense = true;

    run_jacobi_sparse = true;
    run_gauss_seidel_sparse = true;
    run_sor_sparse = true;
    run_cg_sparse = true;
    run_pcg_sparse = true;

    omegas.clear();
    for (int i = 10; i <= 19; ++i) {
        omegas.push_back(static_cast<double>(i) / 10.0);
    }
}

void SolverParam::print_solver_param() const
{
    cout << "================ SolverParam ================" << endl;
    cout << "tol      = " << tol << endl;
    cout << "max_iter = " << max_iter << endl;

    cout << "direct dense full              = " << run_direct_dense_full << endl;
    cout << "direct dense block             = " << run_direct_dense_block << endl;
    cout << "direct sparse block to dense   = " << run_direct_sparse_block_to_dense << endl;
    cout << "cholesky dense                 = " << run_cholesky_dense << endl;
    cout << "ldlt dense                     = " << run_ldlt_dense << endl;

    cout << "jacobi dense                   = " << run_jacobi_dense << endl;
    cout << "gauss seidel dense             = " << run_gauss_seidel_dense << endl;
    cout << "sor dense                      = " << run_sor_dense << endl;
    cout << "cg dense                       = " << run_cg_dense << endl;
    cout << "pcg dense                      = " << run_pcg_dense << endl;

    cout << "jacobi sparse                  = " << run_jacobi_sparse << endl;
    cout << "gauss seidel sparse            = " << run_gauss_seidel_sparse << endl;
    cout << "sor sparse                     = " << run_sor_sparse << endl;
    cout << "cg sparse                      = " << run_cg_sparse << endl;
    cout << "pcg sparse                     = " << run_pcg_sparse << endl;

    cout << "omegas                         = ";
    for (size_t i = 0; i < omegas.size(); ++i) {
        cout << omegas[i];
        if (i + 1 < omegas.size()) cout << ", ";
    }
    cout << endl;
}

// ============================================================
// 6. SpeedParam
// ============================================================

SpeedParam::SpeedParam()
{
    run_mode = BOTH_MODE;
    matrix_mode = BOTH_MATRIX;
    assemble_mode = ASSEMBLE_ALL;
    stress_mode = CENTER_STRESS;

    use_sparse_matrix = true;
    use_fast_mesh = true;
    use_fast_assemble = true;
    use_fast_boundary = true;
    use_fast_post = true;
    use_fast_io = true;

    use_precompute_gauss = true;
    use_precompute_edof = true;
    use_static_local_array = true;
    use_parallel_loop = false;

    compare_base_fast = true;
    compare_solver_result = true;

    timer_detail_level = 2;
}

void SpeedParam::print_speed_param() const
{
    cout << "================ SpeedParam ================" << endl;
    cout << "run_mode              = " << get_run_mode_name(run_mode) << endl;
    cout << "matrix_mode           = " << get_matrix_mode_name(matrix_mode) << endl;
    cout << "assemble_mode         = " << get_assemble_mode_name(assemble_mode) << endl;
    cout << "stress_mode           = " << get_stress_mode_name(stress_mode) << endl;

    cout << "use_sparse_matrix     = " << use_sparse_matrix << endl;
    cout << "use_fast_mesh         = " << use_fast_mesh << endl;
    cout << "use_fast_assemble     = " << use_fast_assemble << endl;
    cout << "use_fast_boundary     = " << use_fast_boundary << endl;
    cout << "use_fast_post         = " << use_fast_post << endl;
    cout << "use_fast_io           = " << use_fast_io << endl;
    cout << "use_precompute_gauss  = " << use_precompute_gauss << endl;
    cout << "use_precompute_edof   = " << use_precompute_edof << endl;
    cout << "use_static_local_array= " << use_static_local_array << endl;
    cout << "use_parallel_loop     = " << use_parallel_loop << endl;
    cout << "compare_base_fast     = " << compare_base_fast << endl;
    cout << "compare_solver_result = " << compare_solver_result << endl;
    cout << "timer_detail_level    = " << timer_detail_level << endl;
}

// ============================================================
// 7. NodeInfo
// ============================================================

NodeInfo::NodeInfo()
{
    id = -1;
    x = 0.0;
    y = 0.0;

    flag_left = 0;
    flag_right = 0;
    flag_top = 0;
    flag_sym = 0;
    flag_hole = 0;
}

// ============================================================
// 8. ElementInfo
// ============================================================

ElementInfo::ElementInfo()
{
    id = -1;
    for (int i = 0; i < 4; ++i) {
        node[i] = -1;
    }
    for (int i = 0; i < 8; ++i) {
        dof[i] = -1;
    }
}

// ============================================================
// 9. BoundaryEdge
// ============================================================

BoundaryEdge::BoundaryEdge()
{
    id = -1;
    n1 = -1;
    n2 = -1;
    marker = 0;
}

// ============================================================
// 10. GaussData
// ============================================================

GaussData::GaussData()
{
    xi = 0.0;
    eta = 0.0;
    weight = 0.0;

    for (int i = 0; i < 4; ++i) {
        N[i] = 0.0;
        dN_dxi[i] = 0.0;
        dN_deta[i] = 0.0;
    }
}

// ============================================================
// 11. SolverInfo
// ============================================================

SolverInfo::SolverInfo()
{
    method = "none";
    matrix_type = "none";

    omega = 0.0;
    converged = false;

    iterations = 0;
    time_seconds = 0.0;
    final_residual = std::numeric_limits<double>::max();

    max_displacement = 0.0;
    n_dof = 0;
    n_solve_dof = 0;
}

// ============================================================
// 12. ElementResult
// ============================================================

ElementResult::ElementResult()
{
    elem_id = -1;

    x_center = 0.0;
    y_center = 0.0;

    epsilon_x = 0.0;
    epsilon_y = 0.0;
    gamma_xy = 0.0;

    sigma_x = 0.0;
    sigma_y = 0.0;
    tau_xy = 0.0;
    sigma_vm = 0.0;
}

// ============================================================
// 13. CaseSummary
// ============================================================

CaseSummary::CaseSummary()
{
    ratio = 0.0;
    a = 0.0;
    b = 0.0;

    max_ux = 0.0;
    max_uy = 0.0;
    max_umag = 0.0;

    max_sigma_x = 0.0;
    max_sigma_y = 0.0;
    max_tau_xy = 0.0;
    max_sigma_vm = 0.0;

    max_sigma_vm_elem = -1;
    max_sigma_vm_x = 0.0;
    max_sigma_vm_y = 0.0;

    Kt_x = 0.0;
    Kt_vm = 0.0;
}

// ============================================================
// 14. CompareBaseFast
// ============================================================

CompareBaseFast::CompareBaseFast()
{
    ratio = 0.0;

    n_node = 0;
    n_elem = 0;

    error_U = 0.0;
    error_stress = 0.0;
    error_Kt_vm = 0.0;

    time_base = 0.0;
    time_fast = 0.0;
    speedup = 0.0;
}

// ============================================================
// 15. CompareSolverResult
// ============================================================

CompareSolverResult::CompareSolverResult()
{
    ratio = 0.0;

    method_ref = "none";
    method_compare = "none";

    error_U = 0.0;
    error_max_disp = 0.0;
    error_residual = 0.0;
}

// ============================================================
// 16. ModelData
// ============================================================

ModelData::ModelData()
{
    clear_data();
}

void ModelData::clear_data()
{
    nodes.clear();
    elements.clear();
    edges.clear();

    K_dense.clear();
    F.clear();
    U.clear();
    U_ref.clear();
    U_base.clear();
    U_fast.clear();

    fixed_dof.clear();
    free_dof.clear();
    is_fixed.clear();

    K_coo.clear_data();

    K_csr.n_row = 0;
    K_csr.n_col = 0;
    K_csr.row_ptr.clear();
    K_csr.col_id.clear();
    K_csr.val.clear();
    K_csr.diag_pos.clear();

    Kff_csr.n_row = 0;
    Kff_csr.n_col = 0;
    Kff_csr.row_ptr.clear();
    Kff_csr.col_id.clear();
    Kff_csr.val.clear();
    Kff_csr.diag_pos.clear();

    Ff.clear();
    Uf.clear();

    gauss_data.clear();
    elem_result.clear();
    elem_result_base.clear();
    elem_result_fast.clear();

    direct_solver_info.clear();
    iterative_solver_info.clear();
    all_solver_info.clear();
    iterative_history.clear();

    summary = CaseSummary();

    n_node = 0;
    n_elem = 0;
    n_dof = 0;

    n_left_fixed_dof = 0;
    n_sym_fixed_dof = 0;
    n_right_load_edge = 0;

    total_load_actual = 0.0;
    total_load_theory = 0.0;

    min_detJ = std::numeric_limits<double>::max();
    max_detJ = -std::numeric_limits<double>::max();

    has_negative_detJ = false;
    solve_success = false;

    sparse_nnz = 0;
    sparse_rate = 0.0;

    warning_info.clear();
}

void ModelData::resize_system(int n_node_in)
{
    n_node = n_node_in;
    n_dof = 2 * n_node;

    K_dense.assign(static_cast<size_t>(n_dof) * static_cast<size_t>(n_dof), 0.0);
    F.assign(n_dof, 0.0);
    U.assign(n_dof, 0.0);
    U_ref.assign(n_dof, 0.0);
    U_base.assign(n_dof, 0.0);
    U_fast.assign(n_dof, 0.0);

    is_fixed.assign(n_dof, 0);
    fixed_dof.clear();
    free_dof.clear();
}

void ModelData::add_warning(const string& msg)
{
    warning_info.push_back(msg);
    cout << "Warning: " << msg << endl;
}

void ModelData::print_model_size() const
{
    cout << "================ ModelData ================" << endl;
    cout << "n_node = " << n_node << endl;
    cout << "n_elem = " << n_elem << endl;
    cout << "n_dof  = " << n_dof << endl;
    cout << "n_edge = " << edges.size() << endl;
}

// ============================================================
// 17. 辅助函数
// ============================================================

string get_run_mode_name(RunMode mode)
{
    switch (mode) {
    case BASE_MODE:
        return "BASE_MODE";
    case FAST_MODE:
        return "FAST_MODE";
    case BOTH_MODE:
        return "BOTH_MODE";
    default:
        return "UNKNOWN_RUN_MODE";
    }
}

string get_matrix_mode_name(MatrixMode mode)
{
    switch (mode) {
    case DENSE_MATRIX:
        return "DENSE_MATRIX";
    case SPARSE_MATRIX:
        return "SPARSE_MATRIX";
    case BOTH_MATRIX:
        return "BOTH_MATRIX";
    default:
        return "UNKNOWN_MATRIX_MODE";
    }
}

string get_assemble_mode_name(AssembleMode mode)
{
    switch (mode) {
    case ASSEMBLE_DENSE_BASE:
        return "ASSEMBLE_DENSE_BASE";
    case ASSEMBLE_DENSE_FAST:
        return "ASSEMBLE_DENSE_FAST";
    case ASSEMBLE_COO_CSR:
        return "ASSEMBLE_COO_CSR";
    case ASSEMBLE_ALL:
        return "ASSEMBLE_ALL";
    default:
        return "UNKNOWN_ASSEMBLE_MODE";
    }
}

string get_stress_mode_name(StressMode mode)
{
    switch (mode) {
    case CENTER_STRESS:
        return "CENTER_STRESS";
    case GAUSS_AVERAGE_STRESS:
        return "GAUSS_AVERAGE_STRESS";
    case BOTH_STRESS:
        return "BOTH_STRESS";
    default:
        return "UNKNOWN_STRESS_MODE";
    }
}

double calculate_umag(double ux, double uy)
{
    return std::sqrt(ux * ux + uy * uy);
}

bool is_close(double a, double b, double eps)
{
    return std::fabs(a - b) <= eps;
}
