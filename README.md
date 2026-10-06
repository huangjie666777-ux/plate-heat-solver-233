# plate_heat233

二维散热板稳态温度场计算库（C++17，P1 线性三角形有限元）：
求解分区材料、多种冷却边界下的 −div(k ∇T) = Q，按单位厚度计算，全部 SI 单位。

## 功能

- 读取 Gmsh 2.2 ASCII 平面网格（z=0），支持二节点线段（类型 1）与三节点三角形（类型 2），
  按**第一物理标签**关联区域材料与边界条件；节点 ID 可不连续，单元绕向任意。
- 网格校验：拒绝重复节点/单元 ID、未知节点引用、退化三角形、非流形边、
  不在外边界的边界线；限制 5000 节点 / 10000 三角形。
- 区域材料：正有限导热率 k [W/(m·K)]，有限体积热源 Q [W/m³]。
- 边界条件（按边界标签）：定温 T、向外热通量 q（向外为正）、对流 h 与环境温度 Ta
  （h 正有限，向外热通量 = h(T−Ta)）；未指定的边界自动绝热。
- 自行装配稀疏对称系统（Eigen3 Sparse）；非零定温采用消元法同步调整右端并保持对称，
  不使用罚数或规则差分。
- 每个连通分量必须含定温或正对流边界，否则报“温度不唯一”。
- 求解器：Eigen SimplicialLDLT；求解失败或最大自由节点残差超过调用者容差则报错，不报成功。
- 输出：按原节点 ID 的温度、逐三角形热流 −k∇T、最大残差、热功率平衡
  （定温边界功率由原系统反力恢复），并导出 ASCII VTK（节点温度 + 单元热流矢量）。

## 构建与测试

依赖：g++ 11.4、Make 4.3、Eigen 3.4.0（已置于 `third_party/eigen`，
亦可用 `make EIGEN=/path/to/eigen` 指定系统安装路径）。

```sh
make            # 构建静态库 build/libplate_heat233.a、示例与自测
make test       # 运行自测
make run-example  # 运行双材料示例，生成 build/dual_material.vtk
```

## 使用示例

```cpp
#include "plate_heat233/plate_heat233.hpp"
using namespace plate_heat233;

Mesh mesh = loadGmsh("plate.msh");

SolveOptions opts;
opts.materials[1] = Material{50.0, 0.0};      // 区域标签 1：k=50, 无热源
opts.materials[2] = Material{5.0, 2.0e4};     // 区域标签 2：k=5, Q=2e4 W/m^3
opts.boundaries[1] = DirichletBC{320.0};      // 边界标签 1：定温 320 K
opts.boundaries[2] = ConvectionBC{25.0, 300.0}; // 边界标签 2：对流 h=25, Ta=300
opts.boundaries[3] = NeumannBC{1000.0};       // 边界标签 3：向外热通量 1000 W/m^2
// 其余边界标签不指定 => 绝热
opts.residualTolerance = 1e-8;

SolveResult res = solve(mesh, opts);          // 失败抛 plate_heat233::Error
writeVtk("out.vtk", mesh, res);
```

完整可运行示例见 `examples/dual_material.cpp`（双材料平板：左侧定温、右侧对流、
上下绝热，右半区含体积热源）。

## 代码结构

| 文件 | 职责 |
| --- | --- |
| `include/plate_heat233/plate_heat233.hpp` | 公共 API 与数据类型 |
| `src/mesh.cpp` | Gmsh 2.2 解析与网格拓扑校验 |
| `src/solve.cpp` | 有限元装配、定温消元、连通性检查、求解、残差与功率平衡 |
| `src/vtk.cpp` | ASCII VTK 导出 |
| `examples/dual_material.cpp` | 双材料调用示例（自带网格生成） |
| `tests/selftest.cpp` | 物理正确性与错误处理自测 |

## 说明

- 所有错误通过异常 `plate_heat233::Error` 报告（继承 `std::runtime_error`）。
- 功率平衡：`source = dirichletOut + neumannOut + convectionOut + imbalance`，
  数值收敛时 `imbalance` 接近零。
- 同节点被两个不同温度的定温边界覆盖时明确失败；相同温度则接受。
