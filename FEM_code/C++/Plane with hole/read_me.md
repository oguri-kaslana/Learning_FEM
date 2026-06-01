# 单头文件版本说明

本版本已将原来的多个 `.h` 文件汇总为一个统一头文件：`fem_all.h`。

使用方式：所有 `.cpp` 文件均统一包含：

```cpp
#include "fem_all.h"
```

原来的各个模块 `.cpp` 文件仍然保留，便于分模块维护和调试。

## 单头文件版本编译命令

```bash
g++ main.cpp fem_class.cpp time_record.cpp sparse_matrix.cpp create_node.cpp create_mesh.cpp create_Jacobi.cpp create_B.cpp create_K.cpp boundary_condition.cpp calculate_result.cpp solve_direct.cpp solve_iterative.cpp solve_algebra.cpp write_result.cpp -o fem_q4_half.exe -O2 -std=c++11
```

如果启用 OpenMP：

```bash
g++ main.cpp fem_class.cpp time_record.cpp sparse_matrix.cpp create_node.cpp create_mesh.cpp create_Jacobi.cpp create_B.cpp create_K.cpp boundary_condition.cpp calculate_result.cpp solve_direct.cpp solve_iterative.cpp solve_algebra.cpp write_result.cpp -o fem_q4_half.exe -O2 -std=c++11 -fopenmp -DUSE_OPENMP
```

---

# FEM_Q4_HalfPlate_Ellipse 使用说明

## 1. 项目功能

本项目用于完成有限元大作业题目：**带中心椭圆孔矩形板的四节点四边形等参元有限元分析**。

程序采用 **C++** 完成有限元主体计算，MATLAB 只负责读取 C++ 输出的 csv 文件并绘制网格、变形图和应力云图。

主要功能包括：

1. 建立带中心椭圆孔矩形板的 **1/2 上半板模型**。
2. 使用 **四节点四边形等参元 Q4** 进行有限元离散。
3. 实现节点、单元、边界边、自由度、单元刚度矩阵、总体刚度矩阵的存储与计算。
4. 输出位移、应变、应力和 von Mises 应力。
5. 计算椭圆孔轴比 b/a = 1, 2, 3, 4, 5 时的孔边应力变化。
6. 实现 dense 和 sparse 两类矩阵存储。
7. 实现直接法、Jacobi、高斯-赛德尔、SOR、CG、PCG 等求解器。
8. 输出每一步耗时、残差、迭代次数、最大位移和求解器对比结果。
9. 通过 MATLAB 绘制各工况应力云图和求解器效率对比图。

---

## 2. 为什么使用 1/2 板模型，而不是 1/4 板模型

完整模型坐标范围为：

```text
x ∈ [0, L]
y ∈ [-H/2, H/2]
```

载荷形式为左端固定、右端受均布拉伸载荷。该问题关于水平中线 y = 0 对称，但关于竖直中线 x = L/2 不对称。

原因是：

- 左端是固定边界；
- 右端是受拉边界；
- 左右边界条件不同；
- 因此不能严格取 1/4 模型。

所以本项目采用上半板模型：

```text
x ∈ [0, L]
y ∈ [0, H/2]
```

---

## 3. 1/2 板模型边界条件

### 3.1 左端固定边界

左端 x = 0：

```text
u = 0
v = 0
```

### 3.2 水平对称边界

水平对称边界 y = 0 被椭圆孔切断，因此分为两段：

```text
左侧：x ∈ [0, cx - a], y = 0
右侧：x ∈ [cx + a, L], y = 0
```

对称边界条件为：

```text
v = 0
u 不约束
```

### 3.3 右端受力边界

右端 x = L 施加沿 +x 方向的均布拉伸载荷 q。

由于采用上半板模型，右端总外力为：

```text
F_total_half = q * (H/2) * t
```

对于右端每条边界边，两个端点各分配一致节点力：

```text
Fx_node = q * le * t / 2
```

其中 le 为边界边长度。

### 3.4 上边界与孔边界

上边界 y = H/2 为自由边界。

椭圆孔边界为自由边界。

---

## 4. 椭圆孔面积保持不变的方法

设基准圆孔半径为 R0，圆孔面积为：

```text
A0 = pi * R0^2
```

椭圆孔面积为：

```text
A = pi * a * b
```

题目要求 b/a 从 1 变化到 5，并保持面积不变。令：

```text
ratio = b/a
```

则：

```text
a = R0 / sqrt(ratio)
b = R0 * sqrt(ratio)
```

这样始终满足：

```text
pi * a * b = pi * R0^2
```

---

## 5. Q4 四节点四边形等参元理论

四节点四边形等参元在自然坐标 xi, eta 中定义。四个形函数为：

```text
N1 = 1/4 * (1 - xi) * (1 - eta)
N2 = 1/4 * (1 + xi) * (1 - eta)
N3 = 1/4 * (1 + xi) * (1 + eta)
N4 = 1/4 * (1 - xi) * (1 + eta)
```

几何映射为：

```text
x = Σ Ni * xi_node
y = Σ Ni * yi_node
```

雅可比矩阵为：

```text
J = [ dx/dxi    dy/dxi
      dx/deta   dy/deta ]
```

通过 J 的逆矩阵可以将自然坐标导数转换为物理坐标导数：

```text
[dN/dx, dN/dy]
```

平面应力问题的 B 矩阵为：

```text
B = [ dN1/dx    0       dN2/dx    0       dN3/dx    0       dN4/dx    0
      0         dN1/dy  0         dN2/dy  0         dN3/dy  0         dN4/dy
      dN1/dy    dN1/dx  dN2/dy    dN2/dx  dN3/dy    dN3/dx  dN4/dy    dN4/dx ]
```

单元刚度矩阵为：

```text
Ke = ∫ B^T D B t dΩ
```

程序使用 2 × 2 高斯积分。

---

## 6. 程序文件结构

```text
FEM_Q4_HalfPlate_Ellipse/
│
├─ main.cpp
│
├─ fem_class.h
├─ fem_class.cpp
│
├─ time_record.h
├─ time_record.cpp
│
├─ sparse_matrix.h
├─ sparse_matrix.cpp
│
├─ create_node.h
├─ create_node.cpp
│
├─ create_mesh.h
├─ create_mesh.cpp
│
├─ create_Jacobi.h
├─ create_Jacobi.cpp
│
├─ create_B.h
├─ create_B.cpp
│
├─ create_K.h
├─ create_K.cpp
│
├─ boundary_condition.h
├─ boundary_condition.cpp
│
├─ calculate_result.h
├─ calculate_result.cpp
│
├─ solve_direct.h
├─ solve_direct.cpp
│
├─ solve_iterative.h
├─ solve_iterative.cpp
│
├─ solve_algebra.h
├─ solve_algebra.cpp
│
├─ write_result.h
├─ write_result.cpp
│
├─ MATLAB_plot/
│  ├─ plot_all.m
│  ├─ plot_mesh.m
│  ├─ plot_deformation.m
│  ├─ plot_stress.m
│  ├─ plot_summary_ratio.m
│  ├─ plot_summary_solver.m
│  └─ plot_summary_time.m
│
├─ output/
│
└─ read_me.md
```

---

## 7. 各文件功能说明

| 文件 | 功能 |
|---|---|
| `main.cpp` | 主程序入口，集中设置参数，循环计算不同 b/a 工况 |
| `fem_class.h/cpp` | 材料、几何、网格、载荷、求解器、节点、单元、结果等数据结构 |
| `time_record.h/cpp` | 计时器类，记录各步骤耗时 |
| `sparse_matrix.h/cpp` | COO、CSR 稀疏矩阵及转换、矩阵向量乘法 |
| `create_node.h/cpp` | 节点创建、边界标记、单元自由度编号 |
| `create_mesh.h/cpp` | 上半板射线法网格生成、边界边生成、网格质量检查 |
| `create_Jacobi.h/cpp` | Q4 形函数、高斯点预计算、雅可比矩阵计算 |
| `create_B.h/cpp` | B 矩阵计算 |
| `create_K.h/cpp` | 单元刚度矩阵和总体刚度矩阵组装 |
| `boundary_condition.h/cpp` | 载荷、固定边界、对称边界、自由自由度集合 |
| `calculate_result.h/cpp` | 应变、应力、von Mises、最大值和 Kt 统计 |
| `solve_direct.h/cpp` | 完整矩阵直接法、自由自由度直接法、Cholesky 接口 |
| `solve_iterative.h/cpp` | Jacobi、GS、SOR、CG、PCG 迭代法 |
| `solve_algebra.h/cpp` | dense/sparse 通用线性代数工具函数 |
| `write_result.h/cpp` | 输出 csv、txt 和汇总文件 |
| `MATLAB_plot/*.m` | 读取 output 结果并绘图 |

---

## 8. 编译方法

### 8.1 普通编译

在项目根目录执行：

```bash
g++ main.cpp fem_class.cpp time_record.cpp sparse_matrix.cpp create_node.cpp create_mesh.cpp create_Jacobi.cpp create_B.cpp create_K.cpp boundary_condition.cpp calculate_result.cpp solve_direct.cpp solve_iterative.cpp solve_algebra.cpp write_result.cpp -o fem_q4_half.exe -O2 -std=c++11
```

### 8.2 启用 OpenMP 的编译方式

如果编译器支持 OpenMP，可以使用：

```bash
g++ main.cpp fem_class.cpp time_record.cpp sparse_matrix.cpp create_node.cpp create_mesh.cpp create_Jacobi.cpp create_B.cpp create_K.cpp boundary_condition.cpp calculate_result.cpp solve_direct.cpp solve_iterative.cpp solve_algebra.cpp write_result.cpp -o fem_q4_half.exe -O2 -std=c++11 -fopenmp -DUSE_OPENMP
```

当前版本中 OpenMP 是可选加速接口，不启用也能正常编译运行。

---

## 9. 运行方法

Linux 或 macOS：

```bash
./fem_q4_half.exe
```

Windows：

```bash
fem_q4_half.exe
```

运行完成后，结果写入：

```text
output/
```

---

## 10. main.cpp 中的重要参数

`main.cpp` 开头集中放置所有可调参数。

常用调参项如下：

| 参数 | 含义 | 默认值 |
|---|---|---|
| `L` | 板长 | 84 mm |
| `H` | 完整板高度 | 84 mm |
| `R0` | 面积等效圆孔半径 | 6 mm |
| `E` | 杨氏模量 | 110000 MPa |
| `nu` | 泊松比 | 0.3 |
| `t` | 厚度 | 1 mm |
| `q` | 右端拉伸载荷 | 1 MPa |
| `n_theta` | 上半椭圆周向分段数 | 96 |
| `n_inner` | 孔边加密层数 | 8 |
| `n_outer` | 外部过渡层数 | 12 |
| `lambda` | 外相似椭圆放大系数 | 1.8 |
| `ratios` | b/a 工况 | 1,2,3,4,5 |

---

## 11. 关于默认的大规模求解器保护

大作业提示词要求所有求解器都要实现并留出开关。本项目已经实现：

- 完整 dense 高斯消元；
- dense 自由自由度子矩阵直接法；
- sparse 子矩阵转 dense 直接法；
- Cholesky；
- Jacobi；
- Gauss-Seidel；
- SOR；
- CG；
- PCG。

但是默认网格 `n_theta = 96, n_inner = 8, n_outer = 12` 时自由度数量较大，完整 dense 直接法、Jacobi、GS、SOR 会非常慢。

因此 `main.cpp` 中设置了：

```cpp
bool auto_skip_expensive_solver = true;
int expensive_dof_limit = 1500;
```

当自由度超过阈值时，程序会自动跳过非常慢的求解器，主要保留 sparse CG 和 sparse PCG。

如果你需要严格对比所有求解器，可以把：

```cpp
auto_skip_expensive_solver = false;
```

或者先把网格调小，例如：

```cpp
n_theta = 24;
n_inner = 4;
n_outer = 6;
```

用于调试和求解器效率对比。

---

## 12. 输出文件说明

每个 ratio 单独输出一个文件夹，例如：

```text
output/ratio_1/
output/ratio_2/
output/ratio_3/
output/ratio_4/
output/ratio_5/
```

每个文件夹内主要文件如下：

| 文件 | 内容 |
|---|---|
| `nodes.csv` | 节点编号、坐标、边界标记 |
| `elements.csv` | 单元编号和四个节点编号 |
| `displacement.csv` | 节点位移 ux、uy、umag |
| `strain.csv` | 单元中心应变 |
| `stress.csv` | 单元中心应力和 von Mises 应力 |
| `direct_solver_summary.csv` | 直接法求解器信息 |
| `iterative_solver_summary.csv` | 迭代法求解器信息 |
| `case_summary.csv` | 单个 ratio 的最大位移、最大应力和 Kt |
| `time_summary.csv` | 单个 ratio 的各步骤耗时 |
| `check_info.txt` | 网格、载荷、detJ、稀疏率和 warning 信息 |

全部工况结束后，在 `output/` 根目录输出：

| 文件 | 内容 |
|---|---|
| `summary_ratio.csv` | 所有 ratio 的最大应力和 Kt 汇总 |
| `summary_solver.csv` | 所有 ratio 的求解器对比 |
| `summary_time.csv` | 所有 ratio 的步骤耗时汇总 |
| `compare_base_fast.csv` | 基础版和加速版耗时对比 |
| `compare_solver_result.csv` | 不同求解器结果指标对比 |

---

## 13. MATLAB 可视化方法

先运行 C++ 程序，生成 `output/` 文件夹。

然后在 MATLAB 中进入：

```matlab
cd MATLAB_plot
plot_all
```

MATLAB 会自动绘制：

### 单个 ratio 图像

- `mesh.png`
- `deformation.png`
- `sigma_x.png`
- `sigma_y.png`
- `tau_xy.png`
- `sigma_vm.png`

### 全部 ratio 汇总图像

- `max_sigma_vm_vs_ratio.png`
- `Kt_vm_vs_ratio.png`
- `max_disp_vs_ratio.png`
- `solver_time_compare.png`
- `iterative_solver_iteration_compare.png`
- `sor_omega_compare.png`
- `cg_pcg_compare.png`
- `dense_sparse_time_compare.png`
- `base_fast_speedup_compare.png`
- `step_time_bar.png`

---

## 14. 基础版本和加速版本的区别

基础版本主要体现最直接、最容易理解的有限元实现流程：

1. 逐节点生成网格；
2. 逐单元计算单刚；
3. dense 总刚矩阵装配；
4. dense 置一法施加边界；
5. 单元中心逐个后处理；
6. 逐行写入结果。

加速版本主要包括：

1. `vector::reserve` 预分配；
2. 节点编号直接由二维索引映射；
3. 高斯点形函数和导数预计算；
4. 单元 edof 预计算；
5. 小矩阵使用固定数组；
6. COO 快速装配；
7. CSR 用于矩阵向量乘法；
8. 自由自由度子矩阵求解；
9. CG / PCG 加速迭代；
10. 快速 I/O；
11. 可选 OpenMP 并行接口。

---

## 15. dense matrix 与 sparse matrix 的区别

### dense matrix

优点：

- 实现简单；
- 便于调试；
- 适合小规模问题。

缺点：

- 内存消耗大；
- 大多数元素为 0，浪费存储；
- 大规模有限元问题中求解速度很慢。

### sparse matrix

优点：

- 只存储非零元；
- 内存占用低；
- 适合有限元总体刚度矩阵；
- 与 CG、PCG 等迭代法配合效率高。

缺点：

- 实现比 dense 复杂；
- 随机修改矩阵不方便；
- 通常需要 COO、CSR 等不同格式配合使用。

---

## 16. 为什么 COO 用于装配，CSR 用于求解

有限元装配阶段需要不断执行：

```text
K[edof[i], edof[j]] += Ke[i,j]
```

这类操作适合 COO 格式，因为 COO 可以直接追加三元组：

```text
row, col, value
```

装配结束后，会有大量重复的 row-col 项，需要合并并排序。转换为 CSR 后，矩阵向量乘法非常高效：

```text
y = K * x
```

CG、PCG、Jacobi、GS、SOR 都需要反复进行矩阵向量运算或按行访问，因此 CSR 更适合求解阶段。

---

## 17. 为什么自由自由度子矩阵求解更快

完整置一法保留所有自由度，包括固定自由度。固定自由度对应的行列会被修改为单位行列。

自由自由度子矩阵方法只求解未知自由度：

```text
Kff * Uf = Ff
```

相比完整矩阵：

1. 求解规模更小；
2. 不需要处理已知位移自由度；
3. 对迭代法更自然；
4. 对 sparse 矩阵更高效。

---

## 18. 应力结果分析建议

### 18.1 b/a = 1

当 b/a = 1 时，椭圆孔退化为圆孔。孔边应力分布相对对称，应力集中区域比较均匀。

### 18.2 b/a 增大

随着 b/a 增大：

```text
a = R0 / sqrt(b/a) 减小
b = R0 * sqrt(b/a) 增大
```

孔洞变得越来越细长。由于孔边曲率分布发生变化，应力流线绕过孔洞时弯曲更剧烈，孔边局部高应力区域会更加集中。

通常可以观察到：

1. 最大应力随 b/a 增大而增大；
2. 高应力区向曲率较大的孔边区域集中；
3. von Mises 应力云图中孔边颜色梯度更加明显；
4. 应力集中系数 Kt_vm 随轴比增大而上升。

---

## 19. 求解器效率分析建议

### 19.1 Jacobi

Jacobi 每次迭代只使用上一轮结果，收敛速度通常最慢。

### 19.2 Gauss-Seidel

Gauss-Seidel 使用当前轮已经更新的新值，通常比 Jacobi 快。

### 19.3 SOR

SOR 在 Gauss-Seidel 基础上加入松弛因子 omega。合适的 omega 可以显著加速，但 omega 过大可能导致振荡或不收敛。

### 19.4 CG

对于施加边界后的线弹性有限元刚度矩阵，通常可以看作对称正定矩阵。CG 对这类问题很适合，通常比 Jacobi、GS、SOR 更高效。

### 19.5 PCG

PCG 在 CG 基础上加入预条件。本项目使用 Jacobi 对角预条件：

```text
M = diag(K)
```

对角预条件实现简单，通常可以改善收敛速度和稳定性。

---

## 20. 如何根据 time_summary.csv 分析程序瓶颈

重点关注以下步骤耗时：

1. `create_K_dense_base`
2. `create_K_dense_fast`
3. `create_K_sparse_coo`
4. `convert_coo_to_csr`
5. `apply_boundary_fast`
6. `solve_*`
7. `calculate_result_fast`
8. `write_result_total`

如果 `create_K_dense_*` 很慢，说明 dense 总刚装配和存储成本较高。

如果 `solve_*` 很慢，说明主要瓶颈在方程组求解阶段。

如果 `write_result_total` 很慢，说明 I/O 输出量较大，可以减少输出频率或只输出必要字段。

---

## 21. 如何根据 compare_base_fast.csv 说明加速效果

`compare_base_fast.csv` 中包含：

```text
ratio,n_node,n_elem,error_U,error_stress,error_Kt_vm,time_base,time_fast,speedup
```

其中：

```text
speedup = time_base / time_fast
```

如果 speedup > 1，说明加速版更快。

如果误差接近 0 或远小于 1e-6，说明加速方案没有破坏计算结果。

---

## 22. 常见问题

### 22.1 程序运行很慢怎么办

建议先调小网格：

```cpp
n_theta = 24;
n_inner = 4;
n_outer = 6;
```

并保留：

```cpp
auto_skip_expensive_solver = true;
```

### 22.2 为什么默认跳过部分求解器

完整 dense 高斯消元、Jacobi、GS、SOR 在默认大网格下会非常慢。程序已经实现这些求解器，但默认用安全开关避免调试时长时间无响应。

### 22.3 MATLAB 没有图像怎么办

请检查：

1. 是否先运行了 C++ 程序；
2. `output/ratio_1` 等文件夹是否存在；
3. `nodes.csv`、`elements.csv`、`stress.csv` 是否存在；
4. MATLAB 当前目录是否为 `MATLAB_plot`。

---

## 23. 作业报告可写的结论模板

本程序采用四节点四边形等参元建立了带中心椭圆孔矩形板的二维平面应力有限元模型。由于原问题左右边界条件不同，仅关于水平中线对称，因此采用上半板 1/2 模型进行计算。通过保持椭圆孔面积不变，改变轴比 b/a，分析孔边应力云图变化。计算结果表明，随着 b/a 增大，椭圆孔变得更加细长，孔边局部曲率变化更加明显，应力流线绕孔弯曲加剧，最大 von Mises 应力和应力集中系数整体呈增大趋势。求解器效率方面，Jacobi 通常最慢，Gauss-Seidel 有一定改善，SOR 对松弛因子敏感，CG 和 PCG 更适合本问题的对称正定刚度矩阵；稀疏矩阵存储和自由自由度缩减可以显著降低计算规模，是有限元程序加速的关键。

