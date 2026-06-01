#ifndef FEM_ALL_H
#define FEM_ALL_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// ============================================================
// fem_all.h
// ------------------------------------------------------------
// 功能：
// 1. 汇总原工程中的全部头文件声明。
// 2. 使用本版本时，所有 .cpp 只需要包含：#include "fem_all.h"
// 3. 原来的 fem_class.h、time_record.h、sparse_matrix.h、create_*.h、
//    boundary_condition.h、calculate_result.h、solve_*.h、write_result.h
//    已经合并到本文件中。
// ============================================================

// ============================================================
// 来自文件：sparse_matrix.h
// ============================================================
// ============================================================
// sparse_matrix.h
// 功能：
// 1. COO 稀疏矩阵格式，用于有限元快速装配。
// 2. CSR 稀疏矩阵格式，用于矩阵向量乘法和迭代求解。
// 3. COO 转 CSR，并自动合并重复 row-col 项。
// 4. 从 CSR 中提取自由自由度子矩阵。
//
// 注意：
// 本文件不依赖 fem_class.h，避免头文件循环包含。
// ============================================================

// ============================================================
// 1. COO 非零元条目
// ============================================================
struct CooEntry {
    int row;
    int col;
    double value;

    CooEntry();
    CooEntry(int row_in, int col_in, double value_in);
};

// ============================================================
// 2. COO 稀疏矩阵
// ============================================================
class SparseCOO {
public:
    int n_row;
    int n_col;
    std::vector<CooEntry> data;

public:
    SparseCOO();
    SparseCOO(int n_row_in, int n_col_in);

    // 设置矩阵尺寸，同时清空已有非零元
    void resize_matrix(int n_row_in, int n_col_in);

    // 添加一个非零元。COO 阶段允许重复 row-col，转换 CSR 时统一合并。
    void add_value(int row, int col, double value);

    // 清空非零元，但保留矩阵尺寸
    void clear_data();

    // 预分配 COO 存储空间，减少 push_back 扩容开销
    void reserve_data(std::size_t n_data);

    // 返回 COO 条目数量，注意不是合并后的 nnz
    int get_nnz_raw() const;
};

// ============================================================
// 3. CSR 稀疏矩阵
// ============================================================
class SparseCSR {
public:
    int n_row;
    int n_col;

    // row_ptr[i] 到 row_ptr[i + 1] - 1 是第 i 行的非零元
    std::vector<int> row_ptr;
    std::vector<int> col_id;
    std::vector<double> val;

    // diag_pos[i] 保存第 i 行对角元在 val 中的位置。
    // 若某行没有对角元，则为 -1。
    std::vector<int> diag_pos;

public:
    SparseCSR();
    SparseCSR(int n_row_in, int n_col_in);

    // 设置矩阵尺寸，同时清空数据
    void resize_matrix(int n_row_in, int n_col_in);

    // 清空所有数据
    void clear_data();

    // CSR 矩阵向量乘法：y = A * x
    void multiply_vector(const std::vector<double>& x,
                         std::vector<double>& y) const;

    // 建立每一行对角元位置，供 GS/SOR/PCG 快速访问
    void create_diag_pos();

    // 获得第 row 行的对角元；若不存在，则返回 0
    double get_diag_value(int row) const;

    // 获得 A(row,col)；主要用于调试和小规模检查，不建议在大循环中频繁调用
    double get_value(int row, int col) const;

    // 返回合并后的非零元数
    int get_nnz() const;

    // 判断矩阵尺寸和 CSR 数组是否基本一致
    bool check_valid(std::string& msg) const;

    // 打印矩阵基本信息
    void print_info(const std::string& name = "SparseCSR") const;
};

// ============================================================
// 4. 稀疏矩阵函数
// ============================================================

// COO 转 CSR：
// 1. 按 row、col 排序；
// 2. 合并重复 row-col 项；
// 3. 建立 row_ptr、col_id、val；
// 4. 建立 diag_pos。
void convert_coo_csr(const SparseCOO& coo, SparseCSR& csr);

// 从 CSR 中提取自由自由度子矩阵 Kff。
// free_dof 保存原矩阵中的自由自由度编号。
void extract_free_csr(const SparseCSR& K,
                      const std::vector<int>& free_dof,
                      SparseCSR& Kff);

// 从完整载荷向量 F 中提取自由自由度载荷 Ff。
void extract_free_vector(const std::vector<double>& F,
                         const std::vector<int>& free_dof,
                         std::vector<double>& Ff);

// 将自由自由度解 Uf 回填到完整位移向量 U。
// 固定自由度默认保持为 0。
void recover_full_vector(const std::vector<double>& Uf,
                         const std::vector<int>& free_dof,
                         int n_dof,
                         std::vector<double>& U);

// 将 CSR 转为 dense。主要用于 sparse block to dense 直接法。
void convert_csr_dense(const SparseCSR& csr,
                       std::vector<double>& A_dense);

// 计算稀疏率 nnz / (n_row * n_col)
double calculate_sparse_rate(const SparseCSR& csr);


// ============================================================
// 来自文件：fem_class.h
// ============================================================
// 注意：ModelData 中需要使用 SparseCOO 和 SparseCSR。
// sparse_matrix.h 后续会单独生成。

// ============================================================
// 1. 运行模式枚举
// ============================================================

// 基础版、加速版、同时对比
enum RunMode {
    BASE_MODE,
    FAST_MODE,
    BOTH_MODE
};

// 矩阵存储模式
enum MatrixMode {
    DENSE_MATRIX,
    SPARSE_MATRIX,
    BOTH_MATRIX
};

// 总刚装配模式
enum AssembleMode {
    ASSEMBLE_DENSE_BASE,
    ASSEMBLE_DENSE_FAST,
    ASSEMBLE_COO_CSR,
    ASSEMBLE_ALL
};

// 应力后处理模式
enum StressMode {
    CENTER_STRESS,
    GAUSS_AVERAGE_STRESS,
    BOTH_STRESS
};

// ============================================================
// 2. 材料参数类
// ============================================================

class MaterialParam {
public:
    double E;          // 杨氏模量，单位 MPa = N/mm^2
    double nu;         // 泊松比
    double t;          // 板厚，单位 mm
    double D[3][3];    // 平面应力本构矩阵

public:
    MaterialParam();

    MaterialParam(double E_in, double nu_in, double t_in);

    // 根据平面应力公式生成本构矩阵
    void build_material();

    // 打印材料参数
    void print_material() const;
};

// ============================================================
// 3. 几何参数类
// ============================================================

class GeometryParam {
public:
    double L;       // 板长，单位 mm
    double H;       // 板高，单位 mm

    double cx;      // 椭圆孔中心 x 坐标
    double cy;      // 椭圆孔中心 y 坐标

    double R0;      // 面积等效圆孔半径
    double ratio;   // 椭圆孔轴比 b/a

    double a;       // 椭圆孔 x 方向半轴
    double b;       // 椭圆孔 y 方向半轴

public:
    GeometryParam();

    GeometryParam(double L_in, double H_in, double R0_in, double ratio_in);

    // 根据 L、H、R0、ratio 计算中心和椭圆半轴
    void build_geometry();

    // 打印几何参数
    void print_geometry() const;
};

// ============================================================
// 4. 网格参数类
// ============================================================

class MeshParam {
public:
    int n_theta;       // 上半椭圆方向分段数，必须为偶数
    int n_inner;       // 椭圆孔到外相似椭圆之间的加密层数
    int n_outer;       // 外相似椭圆到矩形边界的过渡层数

    double lambda;     // 外相似椭圆放大系数

    double min_detJ;   // 网格检查得到的最小 detJ
    double max_detJ;   // 网格检查得到的最大 detJ

    bool use_fast_mesh;

public:
    MeshParam();

    MeshParam(int n_theta_in,
              int n_inner_in,
              int n_outer_in,
              double lambda_in,
              bool use_fast_mesh_in);

    // 打印网格参数
    void print_mesh_param() const;
};

// ============================================================
// 5. 载荷参数类
// ============================================================

class LoadParam {
public:
    double q;             // 均布拉伸载荷，单位 MPa = N/mm^2
    int direction;        // 0 表示 x 方向，1 表示 y 方向
    double total_force;   // 实际装配后的总外力

public:
    LoadParam();

    LoadParam(double q_in, int direction_in);

    // 打印载荷参数
    void print_load() const;
};

// ============================================================
// 6. 求解器参数类
// ============================================================

class SolverParam {
public:
    double tol;       // 迭代收敛容差
    int max_iter;     // 最大迭代次数

    // dense 直接法
    bool run_direct_dense_full;
    bool run_direct_dense_block;
    bool run_direct_sparse_block_to_dense;
    bool run_cholesky_dense;
    bool run_ldlt_dense;

    // dense 迭代法
    bool run_jacobi_dense;
    bool run_gauss_seidel_dense;
    bool run_sor_dense;
    bool run_cg_dense;
    bool run_pcg_dense;

    // sparse 迭代法
    bool run_jacobi_sparse;
    bool run_gauss_seidel_sparse;
    bool run_sor_sparse;
    bool run_cg_sparse;
    bool run_pcg_sparse;

    // SOR 松弛因子列表
    std::vector<double> omegas;

public:
    SolverParam();

    // 打印求解器参数
    void print_solver_param() const;
};

// ============================================================
// 7. 加速参数类
// ============================================================

class SpeedParam {
public:
    RunMode run_mode;
    MatrixMode matrix_mode;
    AssembleMode assemble_mode;
    StressMode stress_mode;

    bool use_sparse_matrix;
    bool use_fast_mesh;
    bool use_fast_assemble;
    bool use_fast_boundary;
    bool use_fast_post;
    bool use_fast_io;

    bool use_precompute_gauss;
    bool use_precompute_edof;
    bool use_static_local_array;
    bool use_parallel_loop;

    bool compare_base_fast;
    bool compare_solver_result;

    int timer_detail_level;

public:
    SpeedParam();

    // 打印加速参数
    void print_speed_param() const;
};

// ============================================================
// 8. 节点信息
// ============================================================

struct NodeInfo {
    int id;       // 节点编号，从 0 开始

    double x;     // 节点 x 坐标
    double y;     // 节点 y 坐标

    int flag_left;    // 左固定边界节点
    int flag_right;   // 右端受力边界节点
    int flag_top;     // 上自由边界节点
    int flag_sym;     // 水平对称边界节点
    int flag_hole;    // 椭圆孔自由边界节点

    NodeInfo();
};

// ============================================================
// 9. 单元信息
// ============================================================

struct ElementInfo {
    int id;       // 单元编号，从 0 开始

    int node[4];  // Q4 单元 4 个节点编号，要求逆时针排列
    int dof[8];   // 单元 8 个自由度编号

    ElementInfo();
};

// ============================================================
// 10. 边界边信息
// ============================================================

struct BoundaryEdge {
    int id;       // 边界边编号
    int n1;       // 边界边起点节点
    int n2;       // 边界边终点节点
    int marker;   // 边界类型标记

    /*
        marker 含义：
        1 = left_fixed
        2 = right_traction
        3 = top_free
        4 = symmetry
        5 = hole_free
    */

    BoundaryEdge();
};

// ============================================================
// 11. 高斯点预计算数据
// ============================================================

struct GaussData {
    double xi;
    double eta;
    double weight;

    double N[4];
    double dN_dxi[4];
    double dN_deta[4];

    GaussData();
};

// ============================================================
// 12. 求解器结果信息
// ============================================================

struct SolverInfo {
    std::string method;       // 求解方法名称
    std::string matrix_type;  // dense 或 sparse

    double omega;             // SOR 松弛因子，其余方法取 0
    bool converged;           // 是否收敛

    int iterations;           // 迭代次数，直接法可取 1
    double time_seconds;      // 求解耗时
    double final_residual;    // 最终相对残差

    double max_displacement;  // 最大位移模长
    int n_dof;                // 总自由度数
    int n_solve_dof;          // 实际求解自由度数

    SolverInfo();
};

// ============================================================
// 12.1 迭代法残差历史
// ============================================================
// 用于补充绘制：残差-迭代步曲线、残差-时间曲线。
struct IterHistory {
    std::string method;
    std::string matrix_type;
    double omega;

    std::vector<int> iteration;
    std::vector<double> residual;
    std::vector<double> time_seconds;

    IterHistory() : method(""), matrix_type(""), omega(0.0) {}
};

// ============================================================
// 13. 单元后处理结果
// ============================================================

struct ElementResult {
    int elem_id;

    double x_center;
    double y_center;

    double epsilon_x;
    double epsilon_y;
    double gamma_xy;

    double sigma_x;
    double sigma_y;
    double tau_xy;
    double sigma_vm;

    ElementResult();
};

// ============================================================
// 14. 单个 ratio 工况汇总结果
// ============================================================

struct CaseSummary {
    double ratio;
    double a;
    double b;

    double max_ux;
    double max_uy;
    double max_umag;

    double max_sigma_x;
    double max_sigma_y;
    double max_tau_xy;
    double max_sigma_vm;

    int max_sigma_vm_elem;
    double max_sigma_vm_x;
    double max_sigma_vm_y;

    double Kt_x;
    double Kt_vm;

    CaseSummary();
};

// ============================================================
// 15. 基础版与加速版对比结果
// ============================================================

struct CompareBaseFast {
    double ratio;

    int n_node;
    int n_elem;

    double error_U;
    double error_stress;
    double error_Kt_vm;

    double time_base;
    double time_fast;
    double speedup;

    CompareBaseFast();
};

// ============================================================
// 16. 不同求解器之间的结果对比
// ============================================================

struct CompareSolverResult {
    double ratio;

    std::string method_ref;
    std::string method_compare;

    double error_U;
    double error_max_disp;
    double error_residual;

    CompareSolverResult();
};

// ============================================================
// 17. 模型数据总类
// ============================================================

class ModelData {
public:
    // ------------------------------
    // 网格数据
    // ------------------------------
    std::vector<NodeInfo> nodes;
    std::vector<ElementInfo> elements;
    std::vector<BoundaryEdge> edges;

    // ------------------------------
    // dense 矩阵与向量
    // ------------------------------
    std::vector<double> K_dense;
    std::vector<double> F;
    std::vector<double> U;

    // 用于保存参考解、基础版结果和加速版结果
    std::vector<double> U_ref;
    std::vector<double> U_base;
    std::vector<double> U_fast;

    // ------------------------------
    // 自由度集合
    // ------------------------------
    std::vector<int> fixed_dof;
    std::vector<int> free_dof;
    std::vector<char> is_fixed;

    // ------------------------------
    // 稀疏矩阵
    // ------------------------------
    SparseCOO K_coo;
    SparseCSR K_csr;

    // 稀疏矩阵自由自由度子矩阵
    SparseCSR Kff_csr;
    std::vector<double> Ff;
    std::vector<double> Uf;

    // ------------------------------
    // 高斯点和结果数据
    // ------------------------------
    std::vector<GaussData> gauss_data;
    std::vector<ElementResult> elem_result;
    std::vector<ElementResult> elem_result_base;
    std::vector<ElementResult> elem_result_fast;

    // ------------------------------
    // 求解器信息
    // ------------------------------
    std::vector<SolverInfo> direct_solver_info;
    std::vector<SolverInfo> iterative_solver_info;
    std::vector<SolverInfo> all_solver_info;

    // 迭代残差历史，用于 MATLAB 绘制 residual-iteration / residual-time 曲线。
    std::vector<IterHistory> iterative_history;

    // ------------------------------
    // 工况汇总
    // ------------------------------
    CaseSummary summary;

    // ------------------------------
    // 检查信息
    // ------------------------------
    int n_node;
    int n_elem;
    int n_dof;

    int n_left_fixed_dof;
    int n_sym_fixed_dof;
    int n_right_load_edge;

    double total_load_actual;
    double total_load_theory;

    double min_detJ;
    double max_detJ;

    bool has_negative_detJ;
    bool solve_success;

    int sparse_nnz;
    double sparse_rate;

    std::vector<std::string> warning_info;

public:
    ModelData();

    // 清空全部数据，进入下一个 ratio 工况前调用
    void clear_data();

    // 根据节点数初始化自由度相关数组
    void resize_system(int n_node_in);

    // 添加 warning 信息
    void add_warning(const std::string& msg);

    // 打印模型基本规模
    void print_model_size() const;
};

// ============================================================
// 18. 辅助函数声明
// ============================================================

// 将 RunMode 转为字符串，方便输出 csv
std::string get_run_mode_name(RunMode mode);

// 将 MatrixMode 转为字符串
std::string get_matrix_mode_name(MatrixMode mode);

// 将 AssembleMode 转为字符串
std::string get_assemble_mode_name(AssembleMode mode);

// 将 StressMode 转为字符串
std::string get_stress_mode_name(StressMode mode);

// 计算节点位移模长
double calculate_umag(double ux, double uy);

// 安全判断浮点数是否接近
bool is_close(double a, double b, double eps = 1.0e-10);


// ============================================================
// 来自文件：time_record.h
// ============================================================
// ============================================================
// TimeRecord
// ------------------------------------------------------------
// 功能：
// 1. 记录有限元程序各主要步骤耗时。
// 2. 支持 start/stop 方式计时。
// 3. 支持手动累加时间。
// 4. 支持输出到屏幕和 csv 文件。
//
// 使用建议：
// 1. 在 main 流程或主要函数入口处 start_time。
// 2. 在主要函数出口处 stop_time。
// 3. 不建议在最底层小函数中频繁计时，避免计时开销影响性能分析。
// ============================================================

class TimeRecord {
public:
    TimeRecord();

    // 开始记录某一步骤时间
    void start_time(const std::string& name);

    // 停止记录某一步骤时间，并自动累加到 time_table 中
    void stop_time(const std::string& name);

    // 手动累加某一步骤时间
    void add_time(const std::string& name, double seconds);

    // 获取某一步骤累计耗时；如果不存在，则返回 0
    double get_time(const std::string& name) const;

    // 判断某一步骤是否已经有计时记录
    bool has_time(const std::string& name) const;

    // 清空全部计时数据
    void clear_time();

    // 打印全部计时结果
    void print_time() const;

    // 输出计时结果到 csv 文件
    // mode 可取：
    // "base"   表示基础版
    // "fast"   表示加速版
    // "total"  表示总流程
    // 其他字符串也会原样写入 mode 列
    void write_time(const std::string& filename, const std::string& mode) const;

    // 获取内部计时表，用于后续汇总所有 ratio 的耗时
    const std::map<std::string, double>& get_table() const;

private:
    std::map<std::string, double> time_table;
    std::map<std::string, std::chrono::high_resolution_clock::time_point> start_table;
};


// ============================================================
// 来自文件：create_node.h
// ============================================================
// ============================================================
// create_node.h
// ------------------------------------------------------------
// 功能：
// 1. 创建单个节点信息。
// 2. 根据坐标和网格索引同步标记边界节点。
// 3. 创建单元自由度编号 edof。
// 4. 统计边界节点数量，便于 check_info.txt 输出。
//
// 说明：
// 真正的网格拓扑在 create_mesh.cpp 中生成，本文件只负责节点
// 信息、边界 flag 和自由度编号等基础操作。
// ============================================================

// 创建单个节点，并根据坐标与索引设置边界标记。
NodeInfo create_node(int node_id,
                     double x,
                     double y,
                     int i_radial,
                     int i_theta,
                     int n_radial,
                     const GeometryParam& geom,
                     const MeshParam& mesh);

// 清空节点全部边界标记。
void clear_node_flag(NodeInfo& node);

// 根据坐标和网格索引设置节点边界标记。
void mark_node(NodeInfo& node,
               int i_radial,
               int i_theta,
               int n_radial,
               const GeometryParam& geom,
               const MeshParam& mesh);

// 仅根据坐标补充节点边界标记。
// 主要用于基础版网格生成或调试场景。
void mark_node_by_coordinate(NodeInfo& node,
                             const GeometryParam& geom,
                             const MeshParam& mesh);

// 生成全部单元的 8 个自由度编号。
void create_edof(ModelData& model,
                 TimeRecord& timer,
                 const std::string& mode_name = "total");

// 更新单个单元的自由度编号。
void update_edof(ElementInfo& elem);

// 节点编号转 u 自由度编号。
int get_dof_u(int node_id);

// 节点编号转 v 自由度编号。
int get_dof_v(int node_id);

// 统计边界节点数量。
void count_boundary_node(const ModelData& model,
                         int& n_left_node,
                         int& n_right_node,
                         int& n_top_node,
                         int& n_sym_node,
                         int& n_hole_node);

// 检查节点编号是否连续、坐标是否有限。
bool check_node(const ModelData& model,
                std::vector<std::string>& warning_info);

// 打印前 max_print 个节点信息，用于调试。
void print_node_info(const ModelData& model,
                     int max_print = 10);

// 判断一个点是否接近某条边界。
bool is_on_left(double x, const GeometryParam& geom, double eps = 1.0e-9);
bool is_on_right(double x, const GeometryParam& geom, double eps = 1.0e-9);
bool is_on_top(double y, const GeometryParam& geom, double eps = 1.0e-9);
bool is_on_bottom(double y, double eps = 1.0e-9);

// 判断节点是否位于水平对称边界的两段有效区间。
bool is_on_symmetry(double x,
                    double y,
                    const GeometryParam& geom,
                    double eps = 1.0e-9);

// 判断节点是否位于椭圆孔边界。
bool is_on_hole(double x,
                double y,
                const GeometryParam& geom,
                double eps = 1.0e-7);


// ============================================================
// 来自文件：create_mesh.h
// ============================================================
// ============================================================
// create_mesh.h
// ------------------------------------------------------------
// 功能：
// 1. 生成带中心椭圆孔矩形板的 1/2 上半板 Q4 网格。
// 2. 采用“内椭圆孔边界 + 外相似椭圆加密带 + 矩形外边界”射线法。
// 3. 同时提供基础版和加速版。
// 4. 生成边界边 edges，用于后续载荷、边界条件和绘图检查。
// 5. 检查 Q4 单元方向和 detJ。
// ============================================================

// 二维点结构，只在网格生成内部使用。
struct MeshPoint {
    double x;
    double y;

    MeshPoint();
    MeshPoint(double x_in, double y_in);
};

// 统一网格生成接口。
void create_mesh(ModelData& model,
                 const GeometryParam& geom,
                 const MeshParam& mesh,
                 const SpeedParam& speed,
                 TimeRecord& timer,
                 const std::string& mode_name = "total");

// 基础版网格生成。
void create_mesh_base(ModelData& model,
                      const GeometryParam& geom,
                      const MeshParam& mesh,
                      TimeRecord& timer,
                      const std::string& mode_name = "base");

// 加速版网格生成。
void create_mesh_fast(ModelData& model,
                      const GeometryParam& geom,
                      const MeshParam& mesh,
                      TimeRecord& timer,
                      const std::string& mode_name = "fast");

// 计算 theta 方向上的内椭圆点。
MeshPoint create_inner_point(double theta,
                             const GeometryParam& geom);

// 计算 theta 方向上的外相似椭圆点。
MeshPoint create_mid_point(double theta,
                           const GeometryParam& geom,
                           const MeshParam& mesh);

// 计算 theta 射线与上半矩形外边界的交点。
MeshPoint create_outer_point(double theta,
                             const GeometryParam& geom,
                             const MeshParam& mesh);

// 在内椭圆、外相似椭圆、矩形外边界之间插值生成径向点。
MeshPoint create_radial_point(int i_radial,
                              int n_radial,
                              const MeshPoint& p_inner,
                              const MeshPoint& p_mid,
                              const MeshPoint& p_outer,
                              const MeshParam& mesh);

// 根据节点网格生成 Q4 单元。
void create_element(ModelData& model,
                    const std::vector<std::vector<int> >& node_grid,
                    const GeometryParam& geom);

// 生成边界边。
void create_boundary_edge(ModelData& model,
                          const std::vector<std::vector<int> >& node_grid,
                          const MeshParam& mesh);

// 检查并修正单元方向。
void check_element_direction(ModelData& model);

// 计算单元中心处的 detJ，用于方向判断和网格质量检查。
double calculate_detJ_center(const ElementInfo& elem,
                             const std::vector<NodeInfo>& nodes);

// 检查网格质量，统计 min_detJ、max_detJ 和负 detJ 警告。
bool check_mesh(ModelData& model,
                MeshParam& mesh,
                TimeRecord& timer,
                const std::string& mode_name = "total");

// 打印网格基本信息。
void print_mesh_info(const ModelData& model);


// ============================================================
// 来自文件：create_Jacobi.h
// ============================================================
// ============================================================
// create_Jacobi.h
// ------------------------------------------------------------
// 功能：
// 1. 创建 Q4 单元 2×2 高斯积分点。
// 2. 计算 Q4 单元雅可比矩阵 J、detJ 和 invJ。
// 3. 提供基础版和加速版接口。
// ============================================================

// 创建 2×2 高斯积分点数据。
void create_Gauss(std::vector<GaussData>& gauss_data,
                  TimeRecord& timer,
                  const std::string& mode_name = "total");

// 根据自然坐标计算 Q4 形函数和自然导数。
void create_shape_Q4(double xi,
                     double eta,
                     double N[4],
                     double dN_dxi[4],
                     double dN_deta[4]);

// 基础版雅可比计算。
bool create_Jacobi_base(const ElementInfo& elem,
                        const std::vector<NodeInfo>& nodes,
                        const GaussData& gp,
                        double J[2][2],
                        double invJ[2][2],
                        double& detJ);

// 加速版雅可比计算。
bool create_Jacobi_fast(const ElementInfo& elem,
                        const std::vector<NodeInfo>& nodes,
                        const GaussData& gp,
                        double J[2][2],
                        double invJ[2][2],
                        double& detJ);

// 统一雅可比接口。
bool create_Jacobi(const ElementInfo& elem,
                   const std::vector<NodeInfo>& nodes,
                   const GaussData& gp,
                   bool use_fast,
                   double J[2][2],
                   double invJ[2][2],
                   double& detJ);

// 在给定自然坐标处计算 detJ，主要用于网格质量检查。
double calculate_detJ(const ElementInfo& elem,
                      const std::vector<NodeInfo>& nodes,
                      double xi,
                      double eta);

// 计算 2×2 矩阵逆矩阵。
bool inverse_matrix_2x2(const double J[2][2],
                        double invJ[2][2],
                        double& detJ);


// ============================================================
// 来自文件：create_B.h
// ============================================================
// ============================================================
// create_B.h
// ------------------------------------------------------------
// 功能：
// 1. 由自然坐标导数和 invJ 计算 dN/dx、dN/dy。
// 2. 构造平面应力 Q4 单元 B 矩阵。
// ============================================================

void create_B_base(const ElementInfo& elem,
                   const std::vector<NodeInfo>& nodes,
                   const GaussData& gp,
                   double B[3][8],
                   double& detJ);

void create_B_fast(const ElementInfo& elem,
                   const std::vector<NodeInfo>& nodes,
                   const GaussData& gp,
                   double B[3][8],
                   double& detJ);

void create_B(const ElementInfo& elem,
              const std::vector<NodeInfo>& nodes,
              const GaussData& gp,
              bool use_fast,
              double B[3][8],
              double& detJ);

void create_B_center(const ElementInfo& elem,
                     const std::vector<NodeInfo>& nodes,
                     bool use_fast,
                     double B[3][8],
                     double& detJ);

void calculate_dN_xy(const GaussData& gp,
                     const double invJ[2][2],
                     double dN_dx[4],
                     double dN_dy[4]);

void clear_B(double B[3][8]);


// ============================================================
// 来自文件：create_K.h
// ============================================================
// ============================================================
// create_K.h
// ------------------------------------------------------------
// 功能：
// 1. 计算 Q4 单元刚度矩阵 Ke。
// 2. 组装 dense 总刚矩阵。
// 3. 组装 COO/CSR 稀疏总刚矩阵。
// 4. 提供基础版、加速版和统一接口。
// ============================================================

void create_K_local(const ElementInfo& elem,
                    const std::vector<NodeInfo>& nodes,
                    const MaterialParam& mat,
                    const std::vector<GaussData>& gauss_data,
                    bool use_fast,
                    double Ke[8][8]);

void create_K_dense_base(ModelData& model,
                         const MaterialParam& mat,
                         TimeRecord& timer,
                         const std::string& mode_name = "base");

void create_K_dense_fast(ModelData& model,
                         const MaterialParam& mat,
                         TimeRecord& timer,
                         const std::string& mode_name = "fast");

void create_K_sparse_coo(ModelData& model,
                         const MaterialParam& mat,
                         TimeRecord& timer,
                         const std::string& mode_name = "sparse");

void create_K_sparse_parallel(ModelData& model,
                              const MaterialParam& mat,
                              const SpeedParam& speed,
                              TimeRecord& timer,
                              const std::string& mode_name = "sparse_parallel");

void create_K_global(ModelData& model,
                     const MaterialParam& mat,
                     const SpeedParam& speed,
                     TimeRecord& timer,
                     const std::string& mode_name = "total");

void clear_local_K(double Ke[8][8]);


// ============================================================
// 来自文件：boundary_condition.h
// ============================================================
// ============================================================
// boundary_condition.h
// ------------------------------------------------------------
// 功能：
// 1. 施加左端固定边界、水平对称边界。
// 2. 施加右端均布拉伸载荷。
// 3. 建立 fixed_dof、free_dof、is_fixed。
// 4. dense 使用置一法，sparse 使用自由自由度子矩阵。
// ============================================================

void apply_fixed_boundary(ModelData& model);
void apply_symmetry_boundary(ModelData& model);

void create_dof_set(ModelData& model,
                    TimeRecord& timer,
                    const std::string& mode_name = "total");

void apply_traction_load(ModelData& model,
                         const GeometryParam& geom,
                         LoadParam& load,
                         TimeRecord& timer,
                         const std::string& mode_name = "total");

void apply_load(ModelData& model,
                const GeometryParam& geom,
                LoadParam& load,
                TimeRecord& timer,
                const std::string& mode_name = "total");

void apply_boundary_dense_base(ModelData& model,
                               TimeRecord& timer,
                               const std::string& mode_name = "base");

void apply_boundary_dense_fast(ModelData& model,
                               TimeRecord& timer,
                               const std::string& mode_name = "fast");

void apply_boundary_sparse(ModelData& model,
                           TimeRecord& timer,
                           const std::string& mode_name = "sparse");

void apply_boundary(ModelData& model,
                    const SpeedParam& speed,
                    TimeRecord& timer,
                    const std::string& mode_name = "total");

bool check_total_load(ModelData& model,
                      const GeometryParam& geom,
                      const MaterialParam& mat,
                      const LoadParam& load,
                      double eps = 1.0e-8);


// ============================================================
// 来自文件：solve_algebra.h
// ============================================================
// ============================================================
// solve_algebra.h
// ------------------------------------------------------------
// dense 和 sparse 求解器共用的基础线性代数函数。
// ============================================================

double dot_vector(const std::vector<double>& a,
                  const std::vector<double>& b);

double norm_vector(const std::vector<double>& a);

void copy_vector(const std::vector<double>& src,
                 std::vector<double>& dst);

void copy_matrix(const std::vector<double>& src,
                 std::vector<double>& dst);

void zero_vector(std::vector<double>& x);

void zero_matrix(std::vector<double>& A);

void multiply_matrix_vector_dense(const std::vector<double>& A,
                                  const std::vector<double>& x,
                                  std::vector<double>& y,
                                  int n);

void multiply_matrix_vector_sparse(const SparseCSR& A,
                                   const std::vector<double>& x,
                                   std::vector<double>& y);

double calculate_residual_dense(const std::vector<double>& A,
                                const std::vector<double>& x,
                                const std::vector<double>& b,
                                int n);

double calculate_residual_sparse(const SparseCSR& A,
                                 const std::vector<double>& x,
                                 const std::vector<double>& b);

void swap_row(std::vector<double>& A,
              std::vector<double>& b,
              int n,
              int r1,
              int r2);

double find_max_abs(const std::vector<double>& x);

double calculate_max_displacement(const std::vector<double>& U);

bool solve_gauss_dense(std::vector<double> A,
                       std::vector<double> b,
                       std::vector<double>& x,
                       int n);

bool solve_cholesky_dense_core(const std::vector<double>& A,
                               const std::vector<double>& b,
                               std::vector<double>& x,
                               int n);

void extract_free_dof_dense(const std::vector<double>& K,
                            const std::vector<double>& F,
                            const std::vector<int>& free_dof,
                            int n_dof,
                            std::vector<double>& Kff,
                            std::vector<double>& Ff);

void extract_free_dof_sparse(const SparseCSR& K,
                             const std::vector<double>& F,
                             const std::vector<int>& free_dof,
                             SparseCSR& Kff,
                             std::vector<double>& Ff);

void convert_sparse_to_dense_submatrix(const SparseCSR& K,
                                       const std::vector<int>& free_dof,
                                       std::vector<double>& Kff_dense);

void extract_diag_dense(const std::vector<double>& A,
                        int n,
                        std::vector<double>& diag);

void extract_diag_sparse(const SparseCSR& A,
                         std::vector<double>& diag);


// ============================================================
// 来自文件：solve_direct.h
// ============================================================
// ============================================================
// solve_direct.h
// ------------------------------------------------------------
// 功能：实现完整 dense、高斯消元、自由自由度子矩阵、sparse 转 dense、Cholesky 等直接法。
// ============================================================

SolverInfo solve_direct_dense_full(ModelData& model,
                                   TimeRecord& timer,
                                   const std::string& mode_name = "dense_full");

SolverInfo solve_direct_dense_block(ModelData& model,
                                    TimeRecord& timer,
                                    const std::string& mode_name = "dense_block");

SolverInfo solve_direct_sparse_block_to_dense(ModelData& model,
                                              TimeRecord& timer,
                                              const std::string& mode_name = "sparse_block_to_dense");

SolverInfo solve_cholesky_dense(ModelData& model,
                                TimeRecord& timer,
                                const std::string& mode_name = "cholesky_dense");

SolverInfo solve_ldlt_dense(ModelData& model,
                            TimeRecord& timer,
                            const std::string& mode_name = "ldlt_dense");

void solve_direct(ModelData& model,
                  const SolverParam& solver,
                  TimeRecord& timer,
                  const std::string& mode_name = "total");


// ============================================================
// 来自文件：solve_iterative.h
// ============================================================
// ============================================================
// solve_iterative.h
// ------------------------------------------------------------
// 功能：实现 Jacobi、Gauss-Seidel、SOR、CG、PCG 的 dense 和 sparse 版本。
// ============================================================

SolverInfo solve_jacobi_dense(ModelData& model, const SolverParam& solver, TimeRecord& timer, const std::string& mode_name = "dense");
SolverInfo solve_gauss_seidel_dense(ModelData& model, const SolverParam& solver, TimeRecord& timer, const std::string& mode_name = "dense");
SolverInfo solve_sor_dense(ModelData& model, const SolverParam& solver, double omega, TimeRecord& timer, const std::string& mode_name = "dense");
SolverInfo solve_cg_dense(ModelData& model, const SolverParam& solver, TimeRecord& timer, const std::string& mode_name = "dense");
SolverInfo solve_pcg_dense(ModelData& model, const SolverParam& solver, TimeRecord& timer, const std::string& mode_name = "dense");

SolverInfo solve_jacobi_sparse(ModelData& model, const SolverParam& solver, TimeRecord& timer, const std::string& mode_name = "sparse");
SolverInfo solve_gauss_seidel_sparse(ModelData& model, const SolverParam& solver, TimeRecord& timer, const std::string& mode_name = "sparse");
SolverInfo solve_sor_sparse(ModelData& model, const SolverParam& solver, double omega, TimeRecord& timer, const std::string& mode_name = "sparse");
SolverInfo solve_cg_sparse(ModelData& model, const SolverParam& solver, TimeRecord& timer, const std::string& mode_name = "sparse");
SolverInfo solve_pcg_sparse(ModelData& model, const SolverParam& solver, TimeRecord& timer, const std::string& mode_name = "sparse");

void solve_iterative(ModelData& model,
                     const SolverParam& solver,
                     TimeRecord& timer,
                     const std::string& mode_name = "total");


// ============================================================
// 来自文件：calculate_result.h
// ============================================================
// ============================================================
// calculate_result.h
// ------------------------------------------------------------
// 功能：计算单元应变、应力、von Mises 应力，并统计最大值和应力集中系数。
// ============================================================

void calculate_strain(const double B[3][8],
                      const double ue[8],
                      double strain[3]);

void calculate_stress(const MaterialParam& mat,
                      const double strain[3],
                      double stress[3]);

double calculate_von_mises(double sigma_x,
                           double sigma_y,
                           double tau_xy);

void calculate_result_base(ModelData& model,
                           const MaterialParam& mat,
                           TimeRecord& timer,
                           const std::string& mode_name = "base");

void calculate_result_fast(ModelData& model,
                           const MaterialParam& mat,
                           TimeRecord& timer,
                           const std::string& mode_name = "fast");

void calculate_result(ModelData& model,
                      const MaterialParam& mat,
                      const LoadParam& load,
                      const SpeedParam& speed,
                      TimeRecord& timer,
                      const std::string& mode_name = "total");

void calculate_result_gauss_average(ModelData& model,
                                    const MaterialParam& mat,
                                    bool use_fast);

void calculate_summary(ModelData& model,
                       const GeometryParam& geom,
                       const LoadParam& load);

void calculate_element_center(const ElementInfo& elem,
                              const std::vector<NodeInfo>& nodes,
                              double& x_center,
                              double& y_center);


// ============================================================
// 来自文件：write_result.h
// ============================================================
// ============================================================
// write_result.h
// ------------------------------------------------------------
// 功能：输出节点、单元、位移、应变、应力、求解器汇总、工况汇总、检查信息和计时信息。
// ============================================================

bool create_folder(const std::string& folder);

void write_nodes(const ModelData& model, const std::string& filename);
void write_elements(const ModelData& model, const std::string& filename);
void write_displacement(const ModelData& model, const std::string& filename);
void write_strain(const ModelData& model, const std::string& filename);
void write_stress(const ModelData& model, const std::string& filename);
void write_case_summary(const ModelData& model, const std::string& filename);
void write_direct_solver_summary(const ModelData& model, const std::string& filename);
void write_iterative_solver_summary(const ModelData& model, const std::string& filename);
void write_check_info(const ModelData& model,
                      const GeometryParam& geom,
                      const MaterialParam& mat,
                      const LoadParam& load,
                      const std::string& filename);

void write_result_base(const ModelData& model,
                       const GeometryParam& geom,
                       const MaterialParam& mat,
                       const LoadParam& load,
                       TimeRecord& timer,
                       const std::string& output_dir);

void write_result_fast(const ModelData& model,
                       const GeometryParam& geom,
                       const MaterialParam& mat,
                       const LoadParam& load,
                       TimeRecord& timer,
                       const std::string& output_dir);

void write_result(const ModelData& model,
                  const GeometryParam& geom,
                  const MaterialParam& mat,
                  const LoadParam& load,
                  const SpeedParam& speed,
                  TimeRecord& timer,
                  const std::string& output_dir);

void write_summary_ratio(const std::vector<CaseSummary>& summaries,
                         const std::string& filename);

void write_summary_solver(const std::vector<SolverInfo>& solver_infos,
                          const std::string& filename);

void write_summary_time(const std::vector<std::pair<std::string, double> >& time_infos,
                        const std::string& filename);

void write_compare_base_fast(const std::vector<CompareBaseFast>& compare_infos,
                             const std::string& filename);

void write_compare_solver_result(const std::vector<CompareSolverResult>& compare_infos,
                                 const std::string& filename);

#endif
