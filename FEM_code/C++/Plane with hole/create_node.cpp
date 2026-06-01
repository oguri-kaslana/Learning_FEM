#include "fem_all.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

using std::cout;
using std::endl;
using std::string;
using std::vector;

// ============================================================
// 1. 创建单个节点
// ============================================================

NodeInfo create_node(int node_id,
                     double x,
                     double y,
                     int i_radial,
                     int i_theta,
                     int n_radial,
                     const GeometryParam& geom,
                     const MeshParam& mesh)
{
    NodeInfo node;
    node.id = node_id;
    node.x = x;
    node.y = y;

    clear_node_flag(node);
    mark_node(node, i_radial, i_theta, n_radial, geom, mesh);

    return node;
}

// ============================================================
// 2. 节点边界标记
// ============================================================

void clear_node_flag(NodeInfo& node)
{
    node.flag_left = 0;
    node.flag_right = 0;
    node.flag_top = 0;
    node.flag_sym = 0;
    node.flag_hole = 0;
}

void mark_node(NodeInfo& node,
               int i_radial,
               int i_theta,
               int n_radial,
               const GeometryParam& geom,
               const MeshParam& mesh)
{
    clear_node_flag(node);

    const double eps = 1.0e-8;

    // i_radial = 0 是上半椭圆孔边界。
    if (i_radial == 0) {
        node.flag_hole = 1;
    }

    // 最外层径向点才可能属于矩形外边界。
    if (i_radial == n_radial - 1) {
        if (is_on_left(node.x, geom, eps)) {
            node.flag_left = 1;
        }
        if (is_on_right(node.x, geom, eps)) {
            node.flag_right = 1;
        }
        if (is_on_top(node.y, geom, eps)) {
            node.flag_top = 1;
        }
    }

    // 水平对称边界只包括 y = 0 的左右两段。
    // 射线网格中 theta = 0 是右侧底边，theta = pi 是左侧底边。
    if ((i_theta == 0 || i_theta == mesh.n_theta) &&
        is_on_symmetry(node.x, node.y, geom, eps)) {
        node.flag_sym = 1;
    }

    // 孔边界自由，不应额外施加对称约束。
    if (node.flag_hole == 1) {
        node.flag_sym = 0;
    }
}

void mark_node_by_coordinate(NodeInfo& node,
                             const GeometryParam& geom,
                             const MeshParam& mesh)
{
    clear_node_flag(node);

    const double eps = 1.0e-8;
    (void)mesh;

    if (is_on_left(node.x, geom, eps)) {
        node.flag_left = 1;
    }
    if (is_on_right(node.x, geom, eps)) {
        node.flag_right = 1;
    }
    if (is_on_top(node.y, geom, eps)) {
        node.flag_top = 1;
    }
    if (is_on_symmetry(node.x, node.y, geom, eps)) {
        node.flag_sym = 1;
    }
    if (is_on_hole(node.x, node.y, geom, 1.0e-6)) {
        node.flag_hole = 1;
        node.flag_sym = 0;
    }
}

// ============================================================
// 3. 自由度编号
// ============================================================

int get_dof_u(int node_id)
{
    return 2 * node_id;
}

int get_dof_v(int node_id)
{
    return 2 * node_id + 1;
}

void update_edof(ElementInfo& elem)
{
    for (int i = 0; i < 4; ++i) {
        const int nid = elem.node[i];
        elem.dof[2 * i] = get_dof_u(nid);
        elem.dof[2 * i + 1] = get_dof_v(nid);
    }
}

void create_edof(ModelData& model,
                 TimeRecord& timer,
                 const string& mode_name)
{
    string step_name = "create_edof";
    if (!mode_name.empty() && mode_name != "total") {
        step_name += "_" + mode_name;
    }

    timer.start_time(step_name);

    for (std::size_t e = 0; e < model.elements.size(); ++e) {
        update_edof(model.elements[e]);
    }

    timer.stop_time(step_name);
}

// ============================================================
// 4. 节点检查和统计
// ============================================================

void count_boundary_node(const ModelData& model,
                         int& n_left_node,
                         int& n_right_node,
                         int& n_top_node,
                         int& n_sym_node,
                         int& n_hole_node)
{
    n_left_node = 0;
    n_right_node = 0;
    n_top_node = 0;
    n_sym_node = 0;
    n_hole_node = 0;

    for (std::size_t i = 0; i < model.nodes.size(); ++i) {
        const NodeInfo& node = model.nodes[i];
        if (node.flag_left) {
            ++n_left_node;
        }
        if (node.flag_right) {
            ++n_right_node;
        }
        if (node.flag_top) {
            ++n_top_node;
        }
        if (node.flag_sym) {
            ++n_sym_node;
        }
        if (node.flag_hole) {
            ++n_hole_node;
        }
    }
}

bool check_node(const ModelData& model,
                vector<string>& warning_info)
{
    bool ok = true;

    for (std::size_t i = 0; i < model.nodes.size(); ++i) {
        const NodeInfo& node = model.nodes[i];

        if (node.id != static_cast<int>(i)) {
            ok = false;
            warning_info.push_back("node id is not continuous at local index " + std::to_string(i));
        }

        if (!std::isfinite(node.x) || !std::isfinite(node.y)) {
            ok = false;
            warning_info.push_back("node coordinate is not finite at node " + std::to_string(node.id));
        }
    }

    return ok;
}

void print_node_info(const ModelData& model,
                     int max_print)
{
    cout << "================ NodeInfo Preview ================" << endl;
    cout << "n_node = " << model.nodes.size() << endl;

    const int n_print = std::min(static_cast<int>(model.nodes.size()), max_print);
    cout << "id, x, y, left, right, top, sym, hole" << endl;

    for (int i = 0; i < n_print; ++i) {
        const NodeInfo& node = model.nodes[i];
        cout << node.id << ", "
             << std::setprecision(8) << node.x << ", "
             << std::setprecision(8) << node.y << ", "
             << node.flag_left << ", "
             << node.flag_right << ", "
             << node.flag_top << ", "
             << node.flag_sym << ", "
             << node.flag_hole << endl;
    }
}

// ============================================================
// 5. 边界判断函数
// ============================================================

bool is_on_left(double x, const GeometryParam& geom, double eps)
{
    return std::fabs(x - 0.0) <= eps * std::max(1.0, geom.L);
}

bool is_on_right(double x, const GeometryParam& geom, double eps)
{
    return std::fabs(x - geom.L) <= eps * std::max(1.0, geom.L);
}

bool is_on_top(double y, const GeometryParam& geom, double eps)
{
    return std::fabs(y - geom.H / 2.0) <= eps * std::max(1.0, geom.H);
}

bool is_on_bottom(double y, double eps)
{
    return std::fabs(y) <= eps;
}

bool is_on_symmetry(double x,
                    double y,
                    const GeometryParam& geom,
                    double eps)
{
    if (!is_on_bottom(y, eps * std::max(1.0, geom.H))) {
        return false;
    }

    // 上半板中，椭圆孔切断了 y=0 对称边界。
    // 有效对称边界为：x <= cx-a 和 x >= cx+a。
    const double left_end = geom.cx - geom.a;
    const double right_start = geom.cx + geom.a;
    const double tol_x = eps * std::max(1.0, geom.L);

    const bool on_left_sym = (x >= -tol_x && x <= left_end + tol_x);
    const bool on_right_sym = (x >= right_start - tol_x && x <= geom.L + tol_x);

    return on_left_sym || on_right_sym;
}

bool is_on_hole(double x,
                double y,
                const GeometryParam& geom,
                double eps)
{
    if (geom.a <= 0.0 || geom.b <= 0.0) {
        return false;
    }

    // 上半椭圆边界：((x-cx)^2/a^2)+((y-cy)^2/b^2)=1, y>=0。
    if (y < -eps) {
        return false;
    }

    const double dx = (x - geom.cx) / geom.a;
    const double dy = (y - geom.cy) / geom.b;
    const double value = dx * dx + dy * dy;

    return std::fabs(value - 1.0) <= eps;
}
