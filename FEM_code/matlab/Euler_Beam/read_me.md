# Euler-Bernoulli 梁有限元程序说明

## 1. 项目简介

本程序用于实现一维 Euler-Bernoulli 梁单元的有限元分析。程序采用二节点 Hermite 梁单元，每个节点包含两个自由度：横向挠度 `w` 和截面转角 `theta`。

程序当前主要用于悬臂梁算例调试：左端固定，右端或全梁施加载荷，并计算节点位移、节点转角、支反力、曲率和弯矩等物理量。

---

## 2. Euler-Bernoulli 梁基本假设

本程序适用于满足 Euler-Bernoulli 梁假设的问题：

1. 梁为细长梁，梁长远大于截面尺寸；
2. 横截面在变形后仍保持平面，并且仍垂直于梁中性轴；
3. 忽略横向剪切变形；
4. 采用小变形、小转角线弹性理论；
5. 材料为线弹性材料，弯曲刚度为 `E*I`。

因此，本程序不适用于厚梁、深梁、剪切变形显著的问题。若需要考虑剪切变形，应改用 Timoshenko 梁单元。

---

## 3. 自由度定义

每个节点有两个自由度：

```matlab
[w, theta]
```

全局位移向量排列为：

```matlab
U_global = [
    w1;
    theta1;
    w2;
    theta2;
    ...
]
```

第 `i` 个节点对应的全局自由度编号为：

```matlab
w_i     -> 2*i - 1
theta_i -> 2*i
```

---

## 4. 文件结构

### 4.1 网格相关函数

| 文件 | 功能 |
|---|---|
| `create_PT.m` | 创建节点坐标矩阵 `P` 和单元拓扑矩阵 `T` |

其中：

```matlab
P(i,:) = [x_i, y_i]
T(e,:) = [node1, node2]
```

---

### 4.2 单元矩阵相关函数

| 文件 | 功能 |
|---|---|
| `calculate_Gauss.m` | 返回 Gauss 积分点和权重 |
| `create_material.m` | 创建弯曲刚度 `E*I` |
| `create_B_local.m` | 创建自然坐标下的 `B_local` |
| `create_Jacobi.m` | 计算一维梁单元雅可比 `J = L_e/2` |
| `create_B_global.m` | 计算每个单元、每个高斯点处的实际坐标 B 矩阵 |
| `create_K_local.m` | 计算所有单元局部刚度矩阵 |
| `create_K_global.m` | 组装整体刚度矩阵 |

欧拉梁单元刚度矩阵采用：

```matlab
Ke = sum(B' * E * I * B * J * weight)
```

当前程序默认使用二点 Gauss 积分。

---

### 4.3 边界条件相关函数

| 文件 | 功能 |
|---|---|
| `boundary_f.m` | 施加节点集中力或集中弯矩 |
| `boundary_q.m` | 使用 Gauss 积分计算均布载荷等效节点力 |
| `boundary_u.m` | 施加位移或转角约束 |

`boundary_f.m` 中的方向编号为：

```matlab
direction = 1  -> w 方向集中力
direction = 2  -> theta 方向集中弯矩
```

`boundary_u.m` 中的方向编号为：

```matlab
direction = 1  -> 约束 w
direction = 2  -> 约束 theta
```

---

### 4.4 求解与后处理函数

| 文件 | 功能 |
|---|---|
| `calculate_u.m` | 求解全局位移向量 |
| `calculate_RF.m` | 计算全局支反力 |
| `calculate_angle.m` | 提取节点转角 |
| `calculate_kappa.m` | 计算各单元各高斯点处的曲率 |
| `calculate_M.m` | 根据曲率计算弯矩 |

支反力计算应使用施加位移边界条件前的原始整体刚度矩阵和原始载荷向量：

```matlab
RF_global = K_global_original * U_global - F_global_original;
```

---

## 5. 主程序说明

主程序文件为：

```matlab
main.m
```

该文件会自动计算三种典型悬臂梁工况：

1. `point_force`：右端集中力；
2. `end_moment`：右端集中弯矩；
3. `uniform_load`：全梁均布载荷。

每个工况都会输出：

1. 节点数、单元数和总自由度数；
2. 右端挠度数值解与解析解；
3. 右端转角数值解与解析解；
4. 左端支反力与支反力矩；
5. 各单元各 Gauss 点处曲率 `kappa`；
6. 各单元各 Gauss 点处弯矩 `M`。

计算完成后，程序会将结果保存为：

```matlab
EulerBeam_debug_results.mat
```

---

## 6. 运行方法

将所有 `.m` 文件放在同一个文件夹下，在 MATLAB 当前路径切换到该文件夹，然后运行：

```matlab
main
```

或直接点击 `main.m` 文件中的运行按钮。

---

## 7. 调试建议

建议优先使用以下参数进行调试：

```matlab
L_total = 1.0;
seeds   = 0.25;
E       = 1.0;
I       = 1.0;
```

这样解析解数值较简单，便于检查程序正确性。程序调通后，再替换为真实材料参数和截面参数。

---

## 8. 注意事项

1. `L_total / seeds` 必须为整数，否则网格划分会出现问题；
2. 当前程序为一维欧拉梁程序，不是二维实体单元程序；
3. `create_Jacobi.m` 中的雅可比是标量 `J = L_e / 2`，不是二维单元的 2×2 雅可比矩阵；
4. `B_global` 是三维数组，其尺寸为：

```matlab
number_elements × number_gauss × 4
```

其中 `B_global(e,g,:)` 表示第 `e` 个单元第 `g` 个 Gauss 点处的 B 矩阵；

5. `boundary_q.m` 适用于均布载荷，若需要线性变化载荷或任意分布载荷，需要修改载荷函数 `q(x)` 的表达方式；
6. 支反力必须使用原始 `K_global` 和原始 `F_global` 计算，不能使用已经施加位移边界条件后的矩阵。

---

## 9. 后续可扩展内容

后续可以在当前程序基础上继续增加：

1. 梁变形曲线绘制函数；
2. 弯矩图绘制函数；
3. 剪力图计算与绘制函数；
4. 多跨梁和不同支座边界条件；
5. 非均布载荷；
6. Timoshenko 梁单元；
7. 稀疏矩阵组装，提高大规模计算效率。
