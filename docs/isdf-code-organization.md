# ISDF 代码分类与变动查看指南

本文是 ISDF 代码的索引，用于快速定位代码和审查变动。

当前代码索引以本文和 `src/GW/ISDF/+isdf/README.md` 为准。
`docs/superpowers/` 下的 plan/spec 是历史开发记录，可能保留旧的
`isdf_*` flat API 名称或迁移前路径。

## 1. 主入口

GW 主流程不直接调用 private helper，而是通过以下层次进入 ISDF：

- `src/GW/epsilon/epsilon_ops.m`：epsilon 的 direct / matrix-elements / reduced-basis 算子选择。
- `src/GW/epsilon/private/*.m`：epsilon full/reduced accumulation 和 page assembly 实现。
- `src/GW/sigma/sigma_ops.m`：sigma 的 matrix element 构造和 screened 收缩选择。
- `src/GW/sigma/private/*.m`：sigma matrix element、full contraction、reduced contraction 实现。
- `src/GW/ISDF/+isdf/*.m`：可复用的 ISDF package API。
- `src/GW/ISDF/+isdf/private/*.m`：只给 package 内部使用的实现细节。

## 2. Public API 分类

Product space:

- `isdf.build_space`
- `isdf.matrix_elements`
- `isdf.real_component`

Reduced polarizability:

- `isdf.polarizability`

Reduced screened interaction:

- `isdf.screened_w`
- `isdf.screened_kernel`
- `isdf.screened_contract`

## 3. Private Helper 分类

Epsilon operators:

- `epsilon_full_init`
- `epsilon_full_evaluate`
- `epsilon_full_accumulate`
- `epsilon_direct_stream_accumulate`
- `epsilon_full_finalize`
- `epsilon_add_state_batch`
- `epsilon_mapped_gme`
- `epsilon_repeat_eden`
- `epsilon_invert_pages`
- `epsilon_reduced_init`
- `epsilon_reduced_evaluate`
- `epsilon_reduced_accumulate`
- `epsilon_reduced_finalize`
- `epsilon_reduced_accumulate_mapped_chi`
- `epsilon_reduced_mapped_zeta_chi`
- `epsilon_page_blkdiag`
- `epsilon_repeat_page_blkdiag`
- `epsilon_gather_screened_w`
- `epsilon_qdata`
- `epsilon_prepared_data`

Sigma operators:

- `sigma_matrix_elements`
- `sigma_contract_full`
- `sigma_contract_reduced`
- `sigma_make_contribution`
- `sigma_kdata`
- `sigma_map_screened_w`
- `sigma_gpu_screened_w`
- `sigma_gather_if_gpu`

Sampling:

- `sample_points`
- `qrcp_sample`
- `scalar_randomized_sample`
- `kmeans_sample`
- `kmeanspp_init`
- `grid_points`
- `distance_to_centers`
- `weighted_sample`
- `weighted_one`
- `unique_fill`

Product-space construction:

- `component_products`
- `component_weight`
- `pair_products`
- `product_gram`
- `zeta_to_g`
- `set_defaults`
- `stable_solve`

Reduced polarizability solver:

- `comega_cstar`
- `comega_direct`
- `comega_cauchy`

Common utility:

- `gather_if_gpu`
- `randn_like`

## 4. 本轮整理后的 Review 路径

如果只想看接口是否变化：

```powershell
git diff -- src/GW/epsilon/epsilon_ops.m src/GW/sigma/sigma_ops.m src/GW/ISDF/+isdf/polarizability.m src/GW/ISDF/+isdf/private/sample_points.m
```

如果想看 epsilon 分类实现：

```powershell
git diff --name-status -- src/GW/epsilon/private
git diff -- src/GW/epsilon/private/epsilon_full*.m src/GW/epsilon/private/epsilon_reduced*.m
```

如果想看 sigma 分类实现：

```powershell
git diff --name-status -- src/GW/sigma/private
git diff -- src/GW/sigma/private/sigma_matrix_elements.m src/GW/sigma/private/sigma_contract*.m
```

如果想看新增 helper 分类：

```powershell
git diff --name-status -- src/GW/ISDF/+isdf/private
```

如果想看 sampling 行为：

```powershell
git diff -- src/GW/ISDF/+isdf/private/sample_points.m src/GW/ISDF/+isdf/private/*sample*.m src/GW/ISDF/+isdf/private/kmeans*.m
```

如果想看 reduced polarizability 行为：

```powershell
git diff -- src/GW/ISDF/+isdf/polarizability.m src/GW/ISDF/+isdf/private/comega*.m
```

如果想看数值防护和 GPU utility：

```powershell
git diff -- src/GW/ISDF/+isdf/private/stable_solve.m src/GW/ISDF/+isdf/private/gather_if_gpu.m src/GW/ISDF/+isdf/private/randn_like.m
```

## 5. 推荐测试

轻量结构和数值测试：

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_comega_cstar.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_component_product_space.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_indices_methods.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_package_api.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_reduced_polarizability.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_screened_kernel_cache.m')"
```

涉及 GPU helper 时再补：

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_gpu_support.m')"
```

## 6. 后续可继续整理的区域

- `epsilon/private`：后续可按 full、reduced、page assembly 再细分 README 或测试。
- `sigma/private`：后续可把 exact static CH 公共逻辑继续提取成 helper。
- `build_space.m`：默认参数已拆到 `private/set_defaults.m`，后续可补更细的 option schema 测试。
- `screened_w.m` 和 `stable_solve.m`：仍有很小的 local cleanup helper，保留是为了让 warning cleanup 靠近调用点。

## 7. 论文中的整体技术实现表述

本节将代码组织映射为可在论文“方法实现”章节中说明的计算流程。它覆盖
GW 主流程、ISDF 降标度和工程实现；ISDF 数学推导及流式内存实现见
`docs/isdf-interpolation-derivation.md`。

### 7.1 端到端计算架构

```text
QE/KSSOLV 数据读取
  → gw_setup：网格、对称性与公共运行时数据
  → epsilon：direct G-space 或 ISDF reduced-basis 屏蔽响应
  → sigma：全空间或 reduced screened-W 自能收缩
  → quasi_energy：准粒子能量
```

- `src/GW/read/` 读入 QE 或 KSSOLV 波函数、能带与相关输入；
- `src/GW/common/runtime/gw_setup.m` 准备公共 GW 数据；
- `epsilon.m` 与 `sigma.m` 分别拥有 q/k/spin 和 band/k/q/spin 的共享
  主循环；
- `epsilon_ops.m` 与 `sigma_ops.m` 在进入循环前选择 direct、
  `matrix_elements` 或 `reduced_basis` 算子，避免在内层循环重复分支；
- `gw_run_json.m` 是可复现生产计算的 JSON 入口。

### 7.2 epsilon 与 sigma 的分层实现

epsilon 的 `epsilon_full_*` 文件构造 conventional G-space 响应，适合作为
物理收敛和 ISDF 误差的参考。`epsilon_reduced_*` 文件在 VC 产品空间中
调用 `isdf.build_space`、`isdf.polarizability` 和 `isdf.screened_w`，形成
约化的 screened interaction。约化极化率通过 `comega_cstar.m` 调度：
`comega_direct.m` 是显式占据--非占据带求和参考，`comega_cauchy.m` 是
静态 Cauchy 降标度求解器。

sigma 的 `sigma_matrix_elements.m` 构造 NN 或 VN 产品空间；
`sigma_contract_full.m` 采用完整介电逆矩阵，
`sigma_contract_reduced.m` 将 screened interaction 投影到目标 ISDF
空间并完成 SEX、COH 与交换相关收缩。NN 是目标带--全部求和带空间，VN
是目标带--占据带空间；后者可在占据带较少时降低裸交换成本。

### 7.3 性能与可扩展性策略

| 策略 | 实现位置 | 作用与代价 |
| --- | --- | --- |
| 随机 QRCP / K-means | `+isdf/private/*sample*.m` | 减少显式带对矩阵构造；随机采样需固定种子并做 rank 收敛。 |
| 实空间波函数缓存 | `real_components.m` | 减少 epsilon/sigma 间重复逆 FFT；增加常驻内存。 |
| global NN/VN 空间 | `sigma_matrix_elements.m` | 多目标带复用采样和辅助基；全局空间更占内存。 |
| screened-kernel 缓存 | `screened_kernel_cache.m` | 避免相同目标空间的重复 W 投影；受字节预算控制。 |
| GPU 数组 | epsilon/sigma/ISDF 核函数 | 利用 GPU 线性代数；GPU 显存也必须纳入 rank 与缓存预算。 |
| 惰性 (P_\mu) 与流式 (C_1) | `product_coefficients.m`、`product_gram.m` | 避免两个大显式矩阵；增加分块循环开销。 |

### 7.4 验证与可复现性

论文应将物理收敛和算法收敛分开报告。前者包括带数、介电 cutoff、q 点
网格和 Coulomb 截断；后者包括 VC/NN/VN rank ratio、采样方法、随机种子、
Cauchy 容差与内存阈值。建议先对可承受的体系使用 `isdf.enable=false` 或
ISDF direct 求解器建立参考，再报告 reduced-basis 结果相对该参考的
epsilon、SX、CH 和准粒子能量误差。

`test/ISDF/unit/` 覆盖产品空间、矩阵元、约化极化率、screened interaction、
GPU 与上下文结构；`test/ISDF/validation/` 提供 AgBr 和 MoS2 的端到端
验证。除物理结果外，论文还应记录 MATLAB 版本、GPU 型号、输入波函数
来源、所有 ISDF 参数、随机种子及缓存/流式开关。
