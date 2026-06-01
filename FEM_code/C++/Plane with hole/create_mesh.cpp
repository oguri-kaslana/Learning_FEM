#include "fem_all.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>

using std::cout;
using std::endl;
using std::string;
using std::vector;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================
// 1. MeshPoint 构造函数
// ============================================================

MeshPoint::MeshPoint()
{
    x = 0.0;
    y = 0.0;
}

MeshPoint::MeshPoint(double x_in, double y_in)
{
    x = x_in;
    y = y_in;
}

// ============================================================
// 2. 统一网格生成接口
// ============================================================

void create_mesh(ModelData& model,
                 const GeometryParam& geom,
                 const MeshParam& mesh,
                 const SpeedParam& speed,
                 TimeRecord& timer,
                 const string& mode_name)
{
    if (speed.use_fast_mesh) {
        create_mesh_fast(model, geom, mesh, timer, mode_name == "total" ? "fast" : mode_name);
    } else {
        create_mesh_base(model, geom, mesh, timer, mode_name == "total" ? "base" : mode_name);
    }
}

// ============================================================
// 3. 基础版网格生成
// ============================================================

void create_mesh_base(ModelData& model,
                      const GeometryParam& geom,
                      const MeshParam& mesh,
                      TimeRecord& timer,
                      const string& mode_name)
{
    string step_name = "create_mesh_base";
    if (!mode_name.empty() && mode_name != "base") {
        step_name += "_" + mode_name;
    }
    timer.start_time(step_name);

    model.nodes.clear();
    model.elements.clear();
    model.edges.clear();

    const int n_theta = mesh.n_theta;
    const int n_radial = mesh.n_inner + mesh.n_outer + 1;

    vector<vector<int> > node_grid(n_radial, vector<int>(n_theta + 1, -1));

    int node_id = 0;

    // 逐点生成节点。基础版不预先 reserve，保留较直观的实现方式。
    for (int i = 0; i < n_radial; ++i) {
        for (int j = 0; j <= n_theta; ++j) {
            const double theta = M_PI * static_cast<double>(j) / static_cast<double>(n_theta);

            MeshPoint p_inner = create_inner_point(theta, geom);
            MeshPoint p_mid = create_mid_point(theta, geom, mesh);
            MeshPoint p_outer = create_outer_point(theta, geom, mesh);
            MeshPoint p = create_radial_point(i, n_radial, p_inner, p_mid, p_outer, mesh);

            NodeInfo node = create_node(node_id, p.x, p.y,
                                        i, j, n_radial,
                                        geom, mesh);

            model.nodes.push_back(node);
            node_grid[i][j] = node_id;
            ++node_id;
        }
    }

    create_element(model, node_grid, geom);
    create_boundary_edge(model, node_grid, mesh);
    check_element_direction(model);

    model.n_node = static_cast<int>(model.nodes.size());
    model.n_elem = static_cast<int>(model.elements.size());
    model.n_dof = 2 * model.n_node;

    timer.stop_time(step_name);
}

// ============================================================
// 4. 加速版网格生成
// ============================================================

void create_mesh_fast(ModelData& model,
                      const GeometryParam& geom,
                      const MeshParam& mesh,
                      TimeRecord& timer,
                      const string& mode_name)
{
    string step_name = "create_mesh_fast";
    if (!mode_name.empty() && mode_name != "fast") {
        step_name += "_" + mode_name;
    }
    timer.start_time(step_name);

    model.nodes.clear();
    model.elements.clear();
    model.edges.clear();

    const int n_theta = mesh.n_theta;
    const int n_radial = mesh.n_inner + mesh.n_outer + 1;
    const int n_node = n_radial * (n_theta + 1);
    const int n_elem = (n_radial - 1) * n_theta;

    model.nodes.reserve(n_node);
    model.elements.reserve(n_elem);
    model.edges.reserve(2 * n_radial + 2 * n_theta + 10);

    vector<vector<int> > node_grid(n_radial, vector<int>(n_theta + 1, -1));

    // 加速点 1：预先计算每条射线的三个控制点。
    vector<MeshPoint> inner_list(n_theta + 1);
    vector<MeshPoint> mid_list(n_theta + 1);
    vector<MeshPoint> outer_list(n_theta + 1);

    for (int j = 0; j <= n_theta; ++j) {
        const double theta = M_PI * static_cast<double>(j) / static_cast<double>(n_theta);
        inner_list[j] = create_inner_point(theta, geom);
        mid_list[j] = create_mid_point(theta, geom, mesh);
        outer_list[j] = create_outer_point(theta, geom, mesh);
    }

    // 加速点 2：节点编号由 i_radial 和 i_theta 直接映射。
    for (int i = 0; i < n_radial; ++i) {
        for (int j = 0; j <= n_theta; ++j) {
            const int node_id = i * (n_theta + 1) + j;
            MeshPoint p = create_radial_point(i, n_radial,
                                              inner_list[j],
                                              mid_list[j],
                                              outer_list[j],
                                              mesh);

            NodeInfo node = create_node(node_id, p.x, p.y,
                                        i, j, n_radial,
                                        geom, mesh);

            model.nodes.push_back(node);
            node_grid[i][j] = node_id;
        }
    }

    create_element(model, node_grid, geom);
    create_boundary_edge(model, node_grid, mesh);
    check_element_direction(model);

    model.n_node = static_cast<int>(model.nodes.size());
    model.n_elem = static_cast<int>(model.elements.size());
    model.n_dof = 2 * model.n_node;

    timer.stop_time(step_name);
}

// ============================================================
// 5. 射线关键点生成
// ============================================================

MeshPoint create_inner_point(double theta,
                             const GeometryParam& geom)
{
    MeshPoint p;
    p.x = geom.cx + geom.a * std::cos(theta);
    p.y = geom.cy + geom.b * std::sin(theta);
    return p;
}

MeshPoint create_mid_point(double theta,
                           const GeometryParam& geom,
                           const MeshParam& mesh)
{
    MeshPoint p;
    p.x = geom.cx + mesh.lambda * geom.a * std::cos(theta);
    p.y = geom.cy + mesh.lambda * geom.b * std::sin(theta);
    return p;
}

MeshPoint create_outer_point(double theta,
                             const GeometryParam& geom,
                             const MeshParam& mesh)
{
    (void)mesh;

    const double eps = 1.0e-12;

    // theta = 0 和 theta = pi 分别对应右、左底边端点。
    if (std::fabs(theta) <= eps) {
        return MeshPoint(geom.L, 0.0);
    }
    if (std::fabs(theta - M_PI) <= eps) {
        return MeshPoint(0.0, 0.0);
    }

    // 射线写成：x = cx + s * a cos(theta), y = cy + s * b sin(theta)。
    // 内椭圆 s = 1，外相似椭圆 s = lambda，外矩形边界 s = s_outer。
    const double dx = geom.a * std::cos(theta);
    const double dy = geom.b * std::sin(theta);

    double best_s = std::numeric_limits<double>::max();

    // 与右边界 x=L 相交。
    if (dx > eps) {
        const double s = (geom.L - geom.cx) / dx;
        if (s > 0.0 && s < best_s) {
            best_s = s;
        }
    }

    // 与左边界 x=0 相交。
    if (dx < -eps) {
        const double s = (0.0 - geom.cx) / dx;
        if (s > 0.0 && s < best_s) {
            best_s = s;
        }
    }

    // 与上边界 y=H/2 相交。
    if (dy > eps) {
        const double s = (geom.H / 2.0 - geom.cy) / dy;
        if (s > 0.0 && s < best_s) {
            best_s = s;
        }
    }

    if (best_s == std::numeric_limits<double>::max()) {
        // 理论上不会出现。若出现，则退化到上边界中点，避免程序崩溃。
        return MeshPoint(geom.cx, geom.H / 2.0);
    }

    MeshPoint p;
    p.x = geom.cx + best_s * dx;
    p.y = geom.cy + best_s * dy;

    // 浮点误差修正，保证外层节点准确落在矩形边界上。
    const double tol = 1.0e-10;
    if (std::fabs(p.x) < tol) {
        p.x = 0.0;
    }
    if (std::fabs(p.x - geom.L) < tol) {
        p.x = geom.L;
    }
    if (std::fabs(p.y) < tol) {
        p.y = 0.0;
    }
    if (std::fabs(p.y - geom.H / 2.0) < tol) {
        p.y = geom.H / 2.0;
    }

    return p;
}

MeshPoint create_radial_point(int i_radial,
                              int n_radial,
                              const MeshPoint& p_inner,
                              const MeshPoint& p_mid,
                              const MeshPoint& p_outer,
                              const MeshParam& mesh)
{
    (void)n_radial;

    MeshPoint p;

    // 0 <= i_radial <= n_inner：椭圆孔到外相似椭圆，加密层。
    if (i_radial <= mesh.n_inner) {
        const double s = static_cast<double>(i_radial) / static_cast<double>(mesh.n_inner);
        p.x = (1.0 - s) * p_inner.x + s * p_mid.x;
        p.y = (1.0 - s) * p_inner.y + s * p_mid.y;
        return p;
    }

    // n_inner < i_radial <= n_inner+n_outer：外相似椭圆到矩形边界。
    const int k = i_radial - mesh.n_inner;
    const double s = static_cast<double>(k) / static_cast<double>(mesh.n_outer);

    p.x = (1.0 - s) * p_mid.x + s * p_outer.x;
    p.y = (1.0 - s) * p_mid.y + s * p_outer.y;

    return p;
}

// ============================================================
// 6. 单元生成与方向检查
// ============================================================

void create_element(ModelData& model,
                    const vector<vector<int> >& node_grid,
                    const GeometryParam& geom)
{
    (void)geom;

    const int n_radial = static_cast<int>(node_grid.size());
    const int n_theta = static_cast<int>(node_grid[0].size()) - 1;

    int elem_id = 0;

    for (int i = 0; i < n_radial - 1; ++i) {
        for (int j = 0; j < n_theta; ++j) {
            ElementInfo elem;
            elem.id = elem_id;

            // 射线网格中的自然排列。
            elem.node[0] = node_grid[i][j];
            elem.node[1] = node_grid[i + 1][j];
            elem.node[2] = node_grid[i + 1][j + 1];
            elem.node[3] = node_grid[i][j + 1];

            update_edof(elem);

            model.elements.push_back(elem);
            ++elem_id;
        }
    }
}

void check_element_direction(ModelData& model)
{
    for (std::size_t e = 0; e < model.elements.size(); ++e) {
        ElementInfo& elem = model.elements[e];
        double detJ = calculate_detJ_center(elem, model.nodes);

        // 若中心点 detJ 为负，交换第 2 和第 4 个节点，使单元方向反转。
        if (detJ < 0.0) {
            std::swap(elem.node[1], elem.node[3]);
            update_edof(elem);
        }
    }
}

double calculate_detJ_center(const ElementInfo& elem,
                             const vector<NodeInfo>& nodes)
{
    // Q4 单元中心 xi=0, eta=0 时的形函数导数。
    const double dN_dxi[4]  = {-0.25,  0.25, 0.25, -0.25};
    const double dN_deta[4] = {-0.25, -0.25, 0.25,  0.25};

    double dx_dxi = 0.0;
    double dy_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_deta = 0.0;

    for (int i = 0; i < 4; ++i) {
        const NodeInfo& node = nodes[elem.node[i]];
        dx_dxi  += dN_dxi[i]  * node.x;
        dy_dxi  += dN_dxi[i]  * node.y;
        dx_deta += dN_deta[i] * node.x;
        dy_deta += dN_deta[i] * node.y;
    }

    return dx_dxi * dy_deta - dy_dxi * dx_deta;
}

// ============================================================
// 7. 边界边生成
// ============================================================

void create_boundary_edge(ModelData& model,
                          const vector<vector<int> >& node_grid,
                          const MeshParam& mesh)
{
    const int n_radial = static_cast<int>(node_grid.size());
    const int n_theta = mesh.n_theta;

    int edge_id = 0;

    // 孔边界：i_radial = 0，marker = 5。
    for (int j = 0; j < n_theta; ++j) {
        BoundaryEdge edge;
        edge.id = edge_id++;
        edge.n1 = node_grid[0][j];
        edge.n2 = node_grid[0][j + 1];
        edge.marker = 5;
        model.edges.push_back(edge);
    }

    // 外边界分段：最外层 i_radial = n_radial - 1。
    // 根据节点 flag 判断是左边、右边还是上边。
    const int i_outer = n_radial - 1;
    for (int j = 0; j < n_theta; ++j) {
        const int n1 = node_grid[i_outer][j];
        const int n2 = node_grid[i_outer][j + 1];
        const NodeInfo& node1 = model.nodes[n1];
        const NodeInfo& node2 = model.nodes[n2];

        BoundaryEdge edge;
        edge.id = edge_id++;
        edge.n1 = n1;
        edge.n2 = n2;
        edge.marker = 3; // 默认上自由边界

        if (node1.flag_right && node2.flag_right) {
            edge.marker = 2;
        } else if (node1.flag_left && node2.flag_left) {
            edge.marker = 1;
        } else if (node1.flag_top && node2.flag_top) {
            edge.marker = 3;
        } else {
            // 射线交到角点附近时，可能形成跨边界小段。该边不施加载荷，作为自由边处理。
            edge.marker = 3;
        }

        model.edges.push_back(edge);
    }

    // 右侧对称边界：theta = 0，从孔边到右端底边，marker = 4。
    for (int i = 0; i < n_radial - 1; ++i) {
        const int n1 = node_grid[i][0];
        const int n2 = node_grid[i + 1][0];

        BoundaryEdge edge;
        edge.id = edge_id++;
        edge.n1 = n1;
        edge.n2 = n2;
        edge.marker = 4;
        model.edges.push_back(edge);
    }

    // 左侧对称边界：theta = pi，从孔边到左端底边，marker = 4。
    for (int i = 0; i < n_radial - 1; ++i) {
        const int n1 = node_grid[i][n_theta];
        const int n2 = node_grid[i + 1][n_theta];

        BoundaryEdge edge;
        edge.id = edge_id++;
        edge.n1 = n1;
        edge.n2 = n2;
        edge.marker = 4;
        model.edges.push_back(edge);
    }
}

// ============================================================
// 8. 网格质量检查
// ============================================================

bool check_mesh(ModelData& model,
                MeshParam& mesh,
                TimeRecord& timer,
                const string& mode_name)
{
    string step_name = "check_mesh";
    if (!mode_name.empty() && mode_name != "total") {
        step_name += "_" + mode_name;
    }
    timer.start_time(step_name);

    bool ok = true;
    double min_detJ = std::numeric_limits<double>::max();
    double max_detJ = -std::numeric_limits<double>::max();
    int n_bad = 0;

    for (std::size_t e = 0; e < model.elements.size(); ++e) {
        const double detJ = calculate_detJ_center(model.elements[e], model.nodes);

        min_detJ = std::min(min_detJ, detJ);
        max_detJ = std::max(max_detJ, detJ);

        if (detJ <= 0.0 || !std::isfinite(detJ)) {
            ok = false;
            ++n_bad;
        }
    }

    mesh.min_detJ = min_detJ;
    mesh.max_detJ = max_detJ;

    model.min_detJ = min_detJ;
    model.max_detJ = max_detJ;
    model.has_negative_detJ = !ok;

    if (!ok) {
        model.add_warning("mesh has element with detJ <= 0, bad element count = " + std::to_string(n_bad));
    }

    int n_left_node = 0;
    int n_right_node = 0;
    int n_top_node = 0;
    int n_sym_node = 0;
    int n_hole_node = 0;
    count_boundary_node(model, n_left_node, n_right_node, n_top_node, n_sym_node, n_hole_node);

    cout << "================ Mesh Check ================" << endl;
    cout << "n_node = " << model.nodes.size() << endl;
    cout << "n_elem = " << model.elements.size() << endl;
    cout << "min_detJ = " << std::setprecision(12) << min_detJ << endl;
    cout << "max_detJ = " << std::setprecision(12) << max_detJ << endl;
    cout << "left nodes = " << n_left_node
         << ", right nodes = " << n_right_node
         << ", top nodes = " << n_top_node
         << ", symmetry nodes = " << n_sym_node
         << ", hole nodes = " << n_hole_node << endl;

    timer.stop_time(step_name);
    return ok;
}

void print_mesh_info(const ModelData& model)
{
    cout << "================ Mesh Info ================" << endl;
    cout << "node number    = " << model.nodes.size() << endl;
    cout << "element number = " << model.elements.size() << endl;
    cout << "edge number    = " << model.edges.size() << endl;
    cout << "dof number     = " << 2 * model.nodes.size() << endl;
}

