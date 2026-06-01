#include "fem_all.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <array>

using std::cout;
using std::endl;
using std::string;
using std::vector;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================
// create_mesh.cpp
// ------------------------------------------------------------
// 本版已由原来的“射线法网格”改为 TFI / Coons patch 网格。
// 思路仿照带孔板结构化网格：将上半板区域分成 4 个 patch，
// 每个 patch 由 4 条边界曲线控制，通过 transfinite interpolation
// 生成协调 Q4 网格。
//
// 当前版本特点：
// 1. 全部单元仍然是标准四节点四边形等参元 Q4；
// 2. 不引入悬挂节点，不引入五节点/六节点过渡单元；
// 3. 径向方向采用 inner_refined 分布：孔边略密，外边界不过度加密；
// 4. 相邻 patch 交界处节点自动合并；
// 5. 通过外边界扫描自动生成左固定、右载荷、上自由、对称、孔边界边。
// ============================================================

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
// 2. 本文件内部辅助函数声明
// ============================================================

namespace {

struct TempElement {
    int node[4];
};

struct PatchKeyPoint {
    MeshPoint P1R, P2R, P3R;
    MeshPoint P4C, P5C;
    MeshPoint CMPR;
    MeshPoint P1L, P2L, P3L;
    MeshPoint CMPL;
};

PatchKeyPoint create_patch_key_point(const GeometryParam& geom)
{
    PatchKeyPoint p;

    const double c45 = std::sqrt(2.0) / 2.0;

    p.P1R = MeshPoint(geom.cx + geom.a, 0.0);
    p.P2R = MeshPoint(geom.L, 0.0);
    p.P3R = MeshPoint(geom.L, geom.H / 2.0);

    p.P4C = MeshPoint(geom.cx, geom.H / 2.0);
    p.P5C = MeshPoint(geom.cx, geom.b);

    p.CMPR = MeshPoint(geom.cx + geom.a * c45, geom.b * c45);

    p.P1L = MeshPoint(geom.cx - geom.a, 0.0);
    p.P2L = MeshPoint(0.0, 0.0);
    p.P3L = MeshPoint(0.0, geom.H / 2.0);

    p.CMPL = MeshPoint(geom.cx - geom.a * c45, geom.b * c45);

    return p;
}

MeshPoint add_point(const MeshPoint& A, const MeshPoint& B)
{
    return MeshPoint(A.x + B.x, A.y + B.y);
}

MeshPoint sub_point(const MeshPoint& A, const MeshPoint& B)
{
    return MeshPoint(A.x - B.x, A.y - B.y);
}

MeshPoint scale_point(double c, const MeshPoint& A)
{
    return MeshPoint(c * A.x, c * A.y);
}

MeshPoint line_point(const MeshPoint& A, const MeshPoint& B, double s)
{
    return MeshPoint((1.0 - s) * A.x + s * B.x,
                     (1.0 - s) * A.y + s * B.y);
}

MeshPoint ellipse_point(const GeometryParam& geom, double theta)
{
    return MeshPoint(geom.cx + geom.a * std::cos(theta),
                     geom.cy + geom.b * std::sin(theta));
}

// 只在孔边加密的径向分布：eta = p^power。
// eta=0 为椭圆孔边界，eta=1 为外边界。
vector<double> create_eta_list(int n_eta_node)
{
    vector<double> eta(n_eta_node, 0.0);

    if (n_eta_node <= 1) {
        eta.assign(2, 0.0);
        eta[1] = 1.0;
        return eta;
    }

    const double power = 1.25;
    for (int i = 0; i < n_eta_node; ++i) {
        const double p = static_cast<double>(i) / static_cast<double>(n_eta_node - 1);
        eta[i] = std::pow(p, power);
    }

    eta.front() = 0.0;
    eta.back() = 1.0;

    return eta;
}

// patch_id：
// 0 = 右下 R1；1 = 右上 R2；2 = 左下 L1；3 = 左上 L2。
MeshPoint get_patch_bottom(int patch_id, double s,
                           const GeometryParam& geom,
                           const PatchKeyPoint& p)
{
    (void)p;

    if (patch_id == 0) {
        return ellipse_point(geom, 0.0 + (M_PI / 4.0) * s);
    }
    if (patch_id == 1) {
        return ellipse_point(geom, M_PI / 2.0 - (M_PI / 4.0) * s);
    }
    if (patch_id == 2) {
        return ellipse_point(geom, M_PI - (M_PI / 4.0) * s);
    }

    return ellipse_point(geom, 3.0 * M_PI / 4.0 - (M_PI / 4.0) * s);
}

MeshPoint get_patch_top(int patch_id, double s,
                        const PatchKeyPoint& p)
{
    if (patch_id == 0) {
        return line_point(p.P2R, p.P3R, s);
    }
    if (patch_id == 1) {
        return line_point(p.P4C, p.P3R, s);
    }
    if (patch_id == 2) {
        return line_point(p.P2L, p.P3L, s);
    }

    return line_point(p.P3L, p.P4C, s);
}

MeshPoint get_patch_left(int patch_id, double t,
                         const PatchKeyPoint& p)
{
    if (patch_id == 0) {
        return line_point(p.P1R, p.P2R, t);
    }
    if (patch_id == 1) {
        return line_point(p.P5C, p.P4C, t);
    }
    if (patch_id == 2) {
        return line_point(p.P1L, p.P2L, t);
    }

    return line_point(p.CMPL, p.P3L, t);
}

MeshPoint get_patch_right(int patch_id, double t,
                          const PatchKeyPoint& p)
{
    if (patch_id == 0) {
        return line_point(p.CMPR, p.P3R, t);
    }
    if (patch_id == 1) {
        return line_point(p.CMPR, p.P3R, t);
    }
    if (patch_id == 2) {
        return line_point(p.CMPL, p.P3L, t);
    }

    return line_point(p.P5C, p.P4C, t);
}

MeshPoint create_coons_point(int patch_id,
                             double s,
                             double t,
                             const GeometryParam& geom,
                             const PatchKeyPoint& key)
{
    const MeshPoint Pb = get_patch_bottom(patch_id, s, geom, key);
    const MeshPoint Pt = get_patch_top(patch_id, s, key);
    const MeshPoint Pl = get_patch_left(patch_id, t, key);
    const MeshPoint Pr = get_patch_right(patch_id, t, key);

    const MeshPoint P00 = get_patch_bottom(patch_id, 0.0, geom, key);
    const MeshPoint P10 = get_patch_bottom(patch_id, 1.0, geom, key);
    const MeshPoint P01 = get_patch_top(patch_id, 0.0, key);
    const MeshPoint P11 = get_patch_top(patch_id, 1.0, key);

    MeshPoint P;

    P.x = (1.0 - t) * Pb.x + t * Pt.x
        + (1.0 - s) * Pl.x + s * Pr.x
        - ((1.0 - s) * (1.0 - t) * P00.x
        + s * (1.0 - t) * P10.x
        + (1.0 - s) * t * P01.x
        + s * t * P11.x);

    P.y = (1.0 - t) * Pb.y + t * Pt.y
        + (1.0 - s) * Pl.y + s * Pr.y
        - ((1.0 - s) * (1.0 - t) * P00.y
        + s * (1.0 - t) * P10.y
        + (1.0 - s) * t * P01.y
        + s * t * P11.y);

    return P;
}

void add_patch(int patch_id,
               int n_seg,
               int n_eta,
               const vector<double>& eta_list,
               const GeometryParam& geom,
               const PatchKeyPoint& key,
               vector<MeshPoint>& temp_nodes,
               vector<TempElement>& temp_elems)
{
    const int nx = n_seg + 1;
    const int ny = n_eta + 1;
    const int node_offset = static_cast<int>(temp_nodes.size());

    vector<vector<int> > grid(ny, vector<int>(nx, -1));

    for (int j = 0; j < ny; ++j) {
        const double t = eta_list[j];
        for (int i = 0; i < nx; ++i) {
            const double s = static_cast<double>(i) / static_cast<double>(n_seg);
            const MeshPoint P = create_coons_point(patch_id, s, t, geom, key);

            const int id = node_offset + static_cast<int>(temp_nodes.size()) - node_offset;
            temp_nodes.push_back(P);
            grid[j][i] = id;
        }
    }

    for (int j = 0; j < ny - 1; ++j) {
        for (int i = 0; i < nx - 1; ++i) {
            TempElement e;
            e.node[0] = grid[j][i];
            e.node[1] = grid[j][i + 1];
            e.node[2] = grid[j + 1][i + 1];
            e.node[3] = grid[j + 1][i];
            temp_elems.push_back(e);
        }
    }
}

double calculate_polygon_area(const vector<NodeInfo>& nodes,
                              const int node_id[4])
{
    double area = 0.0;
    for (int i = 0; i < 4; ++i) {
        const NodeInfo& p1 = nodes[node_id[i]];
        const NodeInfo& p2 = nodes[node_id[(i + 1) % 4]];
        area += p1.x * p2.y - p2.x * p1.y;
    }
    return 0.5 * area;
}

double calculate_polygon_area_temp(const vector<MeshPoint>& nodes,
                                   const int node_id[4])
{
    double area = 0.0;
    for (int i = 0; i < 4; ++i) {
        const MeshPoint& p1 = nodes[node_id[i]];
        const MeshPoint& p2 = nodes[node_id[(i + 1) % 4]];
        area += p1.x * p2.y - p2.x * p1.y;
    }
    return 0.5 * area;
}

bool has_duplicate_node(const int node_id[4])
{
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            if (node_id[i] == node_id[j]) {
                return true;
            }
        }
    }
    return false;
}

void reorder_element_ccw(const vector<NodeInfo>& nodes,
                         int node_id[4])
{
    double cx = 0.0;
    double cy = 0.0;
    for (int i = 0; i < 4; ++i) {
        cx += nodes[node_id[i]].x;
        cy += nodes[node_id[i]].y;
    }
    cx /= 4.0;
    cy /= 4.0;

    vector<std::pair<double, int> > angle_node;
    angle_node.reserve(4);

    for (int i = 0; i < 4; ++i) {
        const NodeInfo& n = nodes[node_id[i]];
        const double angle = std::atan2(n.y - cy, n.x - cx);
        angle_node.push_back(std::make_pair(angle, node_id[i]));
    }

    std::sort(angle_node.begin(), angle_node.end());

    for (int i = 0; i < 4; ++i) {
        node_id[i] = angle_node[i].second;
    }

    if (calculate_polygon_area(nodes, node_id) < 0.0) {
        std::swap(node_id[1], node_id[3]);
    }
}

void merge_patch_nodes(ModelData& model,
                       const vector<MeshPoint>& temp_nodes,
                       const vector<TempElement>& temp_elems,
                       const GeometryParam& geom,
                       const MeshParam& mesh,
                       bool use_reserve)
{
    const double tol = 1.0e-9;

    model.nodes.clear();
    model.elements.clear();
    model.edges.clear();

    if (use_reserve) {
        model.nodes.reserve(temp_nodes.size());
        model.elements.reserve(temp_elems.size());
    }

    std::map<std::pair<long long, long long>, int> node_map;
    vector<int> old_to_new(temp_nodes.size(), -1);

    for (std::size_t i = 0; i < temp_nodes.size(); ++i) {
        const long long kx = static_cast<long long>(std::llround(temp_nodes[i].x / tol));
        const long long ky = static_cast<long long>(std::llround(temp_nodes[i].y / tol));
        const std::pair<long long, long long> key(kx, ky);

        std::map<std::pair<long long, long long>, int>::iterator it = node_map.find(key);

        if (it == node_map.end()) {
            NodeInfo node;
            node.id = static_cast<int>(model.nodes.size());
            node.x = temp_nodes[i].x;
            node.y = temp_nodes[i].y;
            mark_node_by_coordinate(node, geom, mesh);

            const int new_id = node.id;
            node_map[key] = new_id;
            old_to_new[i] = new_id;
            model.nodes.push_back(node);
        } else {
            old_to_new[i] = it->second;
        }
    }

    int elem_id = 0;
    for (std::size_t e = 0; e < temp_elems.size(); ++e) {
        int ids[4];
        for (int i = 0; i < 4; ++i) {
            ids[i] = old_to_new[temp_elems[e].node[i]];
        }

        if (has_duplicate_node(ids)) {
            continue;
        }

        reorder_element_ccw(model.nodes, ids);

        const double area = calculate_polygon_area(model.nodes, ids);
        if (std::fabs(area) < 1.0e-12) {
            continue;
        }

        ElementInfo elem;
        elem.id = elem_id++;
        for (int i = 0; i < 4; ++i) {
            elem.node[i] = ids[i];
        }
        update_edof(elem);
        model.elements.push_back(elem);
    }
}

void create_mesh_tfi(ModelData& model,
                     const GeometryParam& geom,
                     const MeshParam& mesh,
                     bool use_reserve)
{
    int n_theta = mesh.n_theta;
    if (n_theta < 4) {
        n_theta = 4;
    }
    if (n_theta % 4 != 0) {
        const int old_value = n_theta;
        n_theta = ((n_theta + 3) / 4) * 4;
        model.add_warning("TFI mesh requires n_theta divisible by 4; reset n_theta from "
                          + std::to_string(old_value) + " to " + std::to_string(n_theta));
    }

    const int n_seg = n_theta / 4;
    const int n_eta = std::max(2, mesh.n_inner + mesh.n_outer);

    const vector<double> eta_list = create_eta_list(n_eta + 1);
    const PatchKeyPoint key = create_patch_key_point(geom);

    vector<MeshPoint> temp_nodes;
    vector<TempElement> temp_elems;

    if (use_reserve) {
        const int n_patch_node = 4 * (n_seg + 1) * (n_eta + 1);
        const int n_patch_elem = 4 * n_seg * n_eta;
        temp_nodes.reserve(n_patch_node);
        temp_elems.reserve(n_patch_elem);
    }

    // patch 顺序与 MATLAB TFI 脚本保持一致。
    add_patch(0, n_seg, n_eta, eta_list, geom, key, temp_nodes, temp_elems);
    add_patch(1, n_seg, n_eta, eta_list, geom, key, temp_nodes, temp_elems);
    add_patch(2, n_seg, n_eta, eta_list, geom, key, temp_nodes, temp_elems);
    add_patch(3, n_seg, n_eta, eta_list, geom, key, temp_nodes, temp_elems);

    merge_patch_nodes(model, temp_nodes, temp_elems, geom, mesh, use_reserve);
    create_boundary_edge(model, vector<vector<int> >(), mesh);
    check_element_direction(model);

    model.n_node = static_cast<int>(model.nodes.size());
    model.n_elem = static_cast<int>(model.elements.size());
    model.n_dof = 2 * model.n_node;
}

struct EdgeInfo {
    int n1;
    int n2;
    int count;
};

void add_edge_to_map(std::map<std::pair<int, int>, EdgeInfo>& edge_map,
                     int n1,
                     int n2)
{
    const int a = std::min(n1, n2);
    const int b = std::max(n1, n2);
    const std::pair<int, int> key(a, b);

    std::map<std::pair<int, int>, EdgeInfo>::iterator it = edge_map.find(key);
    if (it == edge_map.end()) {
        EdgeInfo info;
        info.n1 = n1;
        info.n2 = n2;
        info.count = 1;
        edge_map[key] = info;
    } else {
        it->second.count += 1;
    }
}

int classify_boundary_edge(const NodeInfo& a,
                           const NodeInfo& b)
{
    if (a.flag_hole && b.flag_hole) {
        return 5;
    }
    if (a.flag_right && b.flag_right) {
        return 2;
    }
    if (a.flag_left && b.flag_left) {
        return 1;
    }
    if (a.flag_sym && b.flag_sym) {
        return 4;
    }
    if (a.flag_top && b.flag_top) {
        return 3;
    }

    // 其他外露边默认视为自由边。
    return 3;
}

double distance_node(const NodeInfo& a,
                     const NodeInfo& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

double calculate_aspect_ratio(const ElementInfo& elem,
                              const vector<NodeInfo>& nodes)
{
    double lmax = 0.0;
    double lmin = std::numeric_limits<double>::max();

    for (int i = 0; i < 4; ++i) {
        const NodeInfo& a = nodes[elem.node[i]];
        const NodeInfo& b = nodes[elem.node[(i + 1) % 4]];
        const double l = distance_node(a, b);
        lmax = std::max(lmax, l);
        lmin = std::min(lmin, l);
    }

    if (lmin <= 1.0e-14) {
        return std::numeric_limits<double>::max();
    }
    return lmax / lmin;
}

void get_shape_derivative_q4(double xi,
                             double eta,
                             double dN_dxi[4],
                             double dN_deta[4])
{
    dN_dxi[0] = -0.25 * (1.0 - eta);
    dN_dxi[1] =  0.25 * (1.0 - eta);
    dN_dxi[2] =  0.25 * (1.0 + eta);
    dN_dxi[3] = -0.25 * (1.0 + eta);

    dN_deta[0] = -0.25 * (1.0 - xi);
    dN_deta[1] = -0.25 * (1.0 + xi);
    dN_deta[2] =  0.25 * (1.0 + xi);
    dN_deta[3] =  0.25 * (1.0 - xi);
}

double calculate_detJ_at(const ElementInfo& elem,
                         const vector<NodeInfo>& nodes,
                         double xi,
                         double eta)
{
    double dN_dxi[4];
    double dN_deta[4];
    get_shape_derivative_q4(xi, eta, dN_dxi, dN_deta);

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

} // namespace

// ============================================================
// 3. 统一网格生成接口
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
// 4. 基础版 TFI 网格生成
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

    // 基础版：不主动 reserve，便于和加速版做时间对比。
    create_mesh_tfi(model, geom, mesh, false);

    timer.stop_time(step_name);
}

// ============================================================
// 5. 加速版 TFI 网格生成
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

    // 加速版：预估节点和单元数量并 reserve，减少 vector 扩容。
    create_mesh_tfi(model, geom, mesh, true);

    timer.stop_time(step_name);
}

// ============================================================
// 6. 兼容旧射线法接口的辅助函数
// ------------------------------------------------------------
// 这些函数保留是为了不破坏 fem_all.h 中已有声明。当前 TFI 网格
// 生成不再依赖它们。
// ============================================================

MeshPoint create_inner_point(double theta,
                             const GeometryParam& geom)
{
    return ellipse_point(geom, theta);
}

MeshPoint create_mid_point(double theta,
                           const GeometryParam& geom,
                           const MeshParam& mesh)
{
    return MeshPoint(geom.cx + mesh.lambda * geom.a * std::cos(theta),
                     geom.cy + mesh.lambda * geom.b * std::sin(theta));
}

MeshPoint create_outer_point(double theta,
                             const GeometryParam& geom,
                             const MeshParam& mesh)
{
    (void)mesh;

    const double eps = 1.0e-12;

    if (std::fabs(theta) <= eps) {
        return MeshPoint(geom.L, 0.0);
    }
    if (std::fabs(theta - M_PI) <= eps) {
        return MeshPoint(0.0, 0.0);
    }

    const double dx = geom.a * std::cos(theta);
    const double dy = geom.b * std::sin(theta);
    double best_s = std::numeric_limits<double>::max();

    if (dx > eps) {
        const double s = (geom.L - geom.cx) / dx;
        const double y = geom.cy + s * dy;
        if (s > 0.0 && y >= -eps && y <= geom.H / 2.0 + eps) {
            best_s = std::min(best_s, s);
        }
    }

    if (dx < -eps) {
        const double s = (0.0 - geom.cx) / dx;
        const double y = geom.cy + s * dy;
        if (s > 0.0 && y >= -eps && y <= geom.H / 2.0 + eps) {
            best_s = std::min(best_s, s);
        }
    }

    if (dy > eps) {
        const double s = (geom.H / 2.0 - geom.cy) / dy;
        const double x = geom.cx + s * dx;
        if (s > 0.0 && x >= -eps && x <= geom.L + eps) {
            best_s = std::min(best_s, s);
        }
    }

    if (best_s == std::numeric_limits<double>::max()) {
        return MeshPoint(geom.cx, geom.H / 2.0);
    }

    MeshPoint p(geom.cx + best_s * dx, geom.cy + best_s * dy);

    const double tol = 1.0e-10;
    if (std::fabs(p.x) < tol) p.x = 0.0;
    if (std::fabs(p.x - geom.L) < tol) p.x = geom.L;
    if (std::fabs(p.y) < tol) p.y = 0.0;
    if (std::fabs(p.y - geom.H / 2.0) < tol) p.y = geom.H / 2.0;

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

    if (i_radial <= mesh.n_inner) {
        const double s = static_cast<double>(i_radial) / static_cast<double>(mesh.n_inner);
        return line_point(p_inner, p_mid, s);
    }

    const int k = i_radial - mesh.n_inner;
    const double s = static_cast<double>(k) / static_cast<double>(mesh.n_outer);
    return line_point(p_mid, p_outer, s);
}

// ============================================================
// 7. 旧结构网格接口保留
// ------------------------------------------------------------
// TFI 网格已经在 create_mesh_tfi 中直接生成单元。该函数保留给
// 旧代码或调试时使用。
// ============================================================

void create_element(ModelData& model,
                    const vector<vector<int> >& node_grid,
                    const GeometryParam& geom)
{
    (void)geom;

    if (node_grid.empty()) {
        return;
    }

    const int n_radial = static_cast<int>(node_grid.size());
    const int n_theta = static_cast<int>(node_grid[0].size()) - 1;

    int elem_id = static_cast<int>(model.elements.size());

    for (int i = 0; i < n_radial - 1; ++i) {
        for (int j = 0; j < n_theta; ++j) {
            ElementInfo elem;
            elem.id = elem_id++;
            elem.node[0] = node_grid[i][j];
            elem.node[1] = node_grid[i + 1][j];
            elem.node[2] = node_grid[i + 1][j + 1];
            elem.node[3] = node_grid[i][j + 1];
            update_edof(elem);
            model.elements.push_back(elem);
        }
    }
}

void check_element_direction(ModelData& model)
{
    for (std::size_t e = 0; e < model.elements.size(); ++e) {
        ElementInfo& elem = model.elements[e];
        double detJ = calculate_detJ_center(elem, model.nodes);

        if (detJ < 0.0) {
            std::swap(elem.node[1], elem.node[3]);
            update_edof(elem);
        }
    }
}

double calculate_detJ_center(const ElementInfo& elem,
                             const vector<NodeInfo>& nodes)
{
    return calculate_detJ_at(elem, nodes, 0.0, 0.0);
}

// ============================================================
// 8. 边界边生成
// ============================================================

void create_boundary_edge(ModelData& model,
                          const vector<vector<int> >& node_grid,
                          const MeshParam& mesh)
{
    (void)node_grid;
    (void)mesh;

    model.edges.clear();

    std::map<std::pair<int, int>, EdgeInfo> edge_map;

    for (std::size_t e = 0; e < model.elements.size(); ++e) {
        const ElementInfo& elem = model.elements[e];
        add_edge_to_map(edge_map, elem.node[0], elem.node[1]);
        add_edge_to_map(edge_map, elem.node[1], elem.node[2]);
        add_edge_to_map(edge_map, elem.node[2], elem.node[3]);
        add_edge_to_map(edge_map, elem.node[3], elem.node[0]);
    }

    int edge_id = 0;
    for (std::map<std::pair<int, int>, EdgeInfo>::const_iterator it = edge_map.begin();
         it != edge_map.end(); ++it) {
        const EdgeInfo& info = it->second;
        if (info.count != 1) {
            continue;
        }

        BoundaryEdge edge;
        edge.id = edge_id++;
        edge.n1 = info.n1;
        edge.n2 = info.n2;
        edge.marker = classify_boundary_edge(model.nodes[edge.n1], model.nodes[edge.n2]);
        model.edges.push_back(edge);
    }
}

// ============================================================
// 9. 网格质量检查
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
    double max_ar = 0.0;
    double sum_ar = 0.0;
    int n_bad_detJ = 0;
    int n_bad_ar = 0;

    const double g = 1.0 / std::sqrt(3.0);
    const double gp[4][2] = {
        {-g, -g}, {g, -g}, {g, g}, {-g, g}
    };

    for (std::size_t e = 0; e < model.elements.size(); ++e) {
        const ElementInfo& elem = model.elements[e];

        for (int k = 0; k < 4; ++k) {
            const double detJ = calculate_detJ_at(elem, model.nodes, gp[k][0], gp[k][1]);
            min_detJ = std::min(min_detJ, detJ);
            max_detJ = std::max(max_detJ, detJ);

            if (detJ <= 0.0 || !std::isfinite(detJ)) {
                ok = false;
                ++n_bad_detJ;
            }
        }

        const double ar = calculate_aspect_ratio(elem, model.nodes);
        max_ar = std::max(max_ar, ar);
        sum_ar += ar;
        if (ar > 8.0) {
            ++n_bad_ar;
        }
    }

    const double mean_ar = model.elements.empty()
                         ? 0.0
                         : sum_ar / static_cast<double>(model.elements.size());

    mesh.min_detJ = min_detJ;
    mesh.max_detJ = max_detJ;

    model.min_detJ = min_detJ;
    model.max_detJ = max_detJ;
    model.has_negative_detJ = !ok;

    if (!ok) {
        model.add_warning("mesh has Gauss point with detJ <= 0, bad Gauss point count = "
                          + std::to_string(n_bad_detJ));
    }

    if (n_bad_ar > 0) {
        model.add_warning("mesh has element aspect ratio > 8, bad element count = "
                          + std::to_string(n_bad_ar)
                          + ", max_aspect_ratio = " + std::to_string(max_ar));
    }

    int n_left_node = 0;
    int n_right_node = 0;
    int n_top_node = 0;
    int n_sym_node = 0;
    int n_hole_node = 0;
    count_boundary_node(model, n_left_node, n_right_node, n_top_node, n_sym_node, n_hole_node);

    cout << "================ Mesh Check ================" << endl;
    cout << "mesh type = TFI / Coons patch, eta mode = inner_refined" << endl;
    cout << "n_node = " << model.nodes.size() << endl;
    cout << "n_elem = " << model.elements.size() << endl;
    cout << "min_detJ = " << std::setprecision(12) << min_detJ << endl;
    cout << "max_detJ = " << std::setprecision(12) << max_detJ << endl;
    cout << "max_aspect_ratio = " << std::setprecision(8) << max_ar << endl;
    cout << "mean_aspect_ratio = " << std::setprecision(8) << mean_ar << endl;
    cout << "bad_aspect_count(AR>8) = " << n_bad_ar << endl;
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
    cout << "n_node = " << model.nodes.size() << endl;
    cout << "n_elem = " << model.elements.size() << endl;
    cout << "n_edge = " << model.edges.size() << endl;
    cout << "n_dof  = " << 2 * model.nodes.size() << endl;
}
