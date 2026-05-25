# Timoshenko 梁有限元程序说明

## 1. 项目简介

本程序用于实现一维二节点 Timoshenko 梁单元的有限元分析。程序采用每个节点两个自由度的格式：横向挠度 `w` 和独立截面转角 `theta`。

程序当前主要用于悬臂梁调试算例：左端固定，右端或全梁施加载荷，并计算节点位移、节点转角、支反力、曲率、弯矩、剪切应变和剪力等物理量。

---

## 2. Timoshenko 梁基本假设

Timoshenko 梁理论适用于需要考虑横向剪切变形的梁问题，尤其适合中等厚度梁或短梁。与 Euler-Bernoulli 梁不同，Timoshenko 梁不再要求截面转角 `theta` 等于挠度斜率 `dw/dx`。

本程序采用以下关系：

```matlab
kappa = dtheta/dx
gamma = dw/dx - theta
M = E * I * kappa
Q = kappa_s * G * A * gamma
```

其中：

| 符号 | 含义 |
|---|---|
| `kappa` | 曲率 |
| `gamma` | 剪切应变 |
| `M` | 弯矩 |
| `Q` | 剪力 |
| `E` | 弹性模量 |
| `G` | 剪切模量 |
| `A` | 截面面积 |
| `I` | 截面惯性矩 |
| `kappa_s` | 剪切修正系数，矩形截面常取 `5/6` |

---

## 3. 与 Euler-Bernoulli 梁程序的主要区别

| 内容 | Euler-Bernoulli 梁 | Timoshenko 梁 |
|---|---|---|
| 截面转角 | `theta = dw/dx` | `theta` 是独立自由度 |
| 曲率 | `kappa = d²w/dx²` | `kappa = dtheta/dx` |
| 剪切变形 | 忽略 | 考虑 |
| 剪切应变 | 不计算 | `gamma = dw/dx - theta` |
| 材料刚度 | `E*I` | `diag(E*I, kappa_s*G*A)` |
| 刚度积分 | 只需弯曲项 | 弯曲项 + 剪切项 |
| 均布载荷形函数 | Hermite 形函数 | 线性位移形函数 `[N1; 0; N2; 0]` |

为减轻 Timoshenko 梁在细长梁情况下可能出现的剪切锁死，本程序采用选择性积分：

```matlab
弯曲刚度项：2 点 Gauss 积分
剪切刚度项：1 点 Gauss 减缩积分
```

---

## 4. 自由度定义

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

## 5. 文件结构

### 5.1 网格函数

| 文件 | 功能 |
|---|---|
| `create_PT.m` | 创建节点坐标矩阵 `P` 和单元拓扑矩阵 `T` |

其中：

```matlab
P(i,:) = [x_i, y_i]
T(e,:) = [node1, node2]
```

---

### 5.2 单元矩阵函数

| 文件 | 功能 |
|---|---|
| `calculate_Gauss.m` | 返回一维 Gauss 积分点和权重 |
| `create_material.m` | 创建 Timoshenko 梁材料刚度矩阵 `diag(EI, kGA)` |
| `create_B_local.m` | 创建弯曲 B 和剪切 B |
| `create_Jacobi.m` | 计算一维梁单元雅可比 `J = L_e/2` |
| `create_B_global.m` | 计算每个单元各 Gauss 点处的 `Bb` 和 `Bs` |
| `create_K_local.m` | 计算所有单元局部刚度矩阵 |
| `create_K_global.m` | 组装整体刚度矩阵 |

单元刚度矩阵由两部分组成：

```matlab
Ke = Ke_bending + Ke_shear
```

其中：

```matlab
Ke_bending = integral(Bb' * E*I * Bb dx)
Ke_shear   = integral(Bs' * kappa_s*G*A * Bs dx)
```

---

### 5.3 边界条件函数

| 文件 | 功能 |
|---|---|
| `boundary_f.m` | 施加节点集中力或集中弯矩 |
| `boundary_q.m` | 使用 Gauss 积分计算均布载荷等效节点力 |
| `boundary_u.m` | 施加位移或转角约束 |

方向编号为：

```matlab
direction = 1  -> w 自由度
direction = 2  -> theta 自由度
```

---

### 5.4 求解与后处理函数

| 文件 | 功能 |
|---|---|
| `calculate_u.m` | 求解全局位移向量 |
| `calculate_RF.m` | 计算全局支反力 |
| `calculate_angle.m` | 提取节点转角 |
| `calculate_kappa.m` | 计算曲率 `kappa = dtheta/dx` |
| `calculate_M.m` | 计算弯矩 `M = E*I*kappa` |
| `calculate_gamma.m` | 计算剪切应变 `gamma = dw/dx - theta` |
| `calculate_Q.m` | 计算剪力 `Q = kappa_s*G*A*gamma` |

支反力计算应使用施加位移边界条件前的原始刚度矩阵和原始载荷向量：

```matlab
RF_global = K_global_original * U_global - F_global_original;
```

---

## 6. 主程序说明

主程序文件为：

```matlab
main.m
```

该文件自动计算三种典型悬臂梁工况：

1. `point_force`：右端集中力；
2. `end_moment`：右端集中弯矩；
3. `uniform_load`：全梁均布载荷。

每个工况会输出：

1. 节点数、单元数和总自由度数；
2. 右端挠度数值解与 Timoshenko 梁解析解；
3. 右端转角数值解与解析解；
4. 左端支反力与支反力矩；
5. 各单元 Gauss 点处曲率 `kappa`；
6. 各单元 Gauss 点处弯矩 `M`；
7. 各单元剪切应变 `gamma`；
8. 各单元剪力 `Q`。

计算完成后，程序会保存：

```matlab
TimoshenkoBeam_debug_results.mat
```

---

## 7. 运行方法

将所有 `.m` 文件放在同一个文件夹下，在 MATLAB 当前路径切换到该文件夹，然后运行：

```matlab
main
```

---

## 8. 注意事项

1. `L_total / seeds` 必须为整数；
2. 当前程序是一维二节点 Timoshenko 梁程序，不是二维实体单元程序；
3. `create_Jacobi.m` 中的雅可比是标量 `J = L_e / 2`；
4. `create_B_global.m` 会分别输出弯曲 B 和剪切 B：

```matlab
[Bb_global, Bs_global] = create_B_global(P)
```

5. 剪切刚度项采用 1 点 Gauss 减缩积分，主要用于减轻剪切锁死；
6. 如果研究非常细长的梁，Timoshenko 梁结果应逐渐接近 Euler-Bernoulli 梁结果；
7. 如果研究较厚梁或短梁，剪切变形对挠度的贡献会更明显。

---

## 9. 后续可扩展内容

后续可以继续增加：

1. 梁变形曲线绘制；
2. 弯矩图和剪力图绘制；
3. 非均布载荷；
4. 多跨梁和弹性支座；
5. 与 Euler-Bernoulli 梁单元的收敛对比；
6. 剪切锁死与减缩积分的数值实验。
