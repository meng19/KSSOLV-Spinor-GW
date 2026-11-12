# KSSOLV Spinor GW Python 迁移方案

## 1. 范围与迁移原则

- MATLAB 仓库 `C:\Users\USTC\Documents\GitHub\KSSOLV-Spinor-GW` 是唯一数值与行为参照，目标仓库为 `F:\KSSOLV_Spinor_GW_python`。
- 本期覆盖 QE/KSSOLV 输入、公共数值核、direct/static COHSEX/full-frequency、scalar/spinor、单 k/多 k、ISDF `matrix_elements`/`reduced_basis`、exact static CH、准粒子能量、检查点及 JSON 工作流。
- 不迁移 SCF、EigSolver、GeomOpt、Tools，不新增 BSE 功能；仅保留 GW 所需 kernel helper。
- 保留目标仓库现有包名和可复用代码，不做一次性全目录重写。每一阶段通过回归后再进入下一阶段。
- CPU 首先完整落地：NumPy 负责数组/FFT，SciPy 负责 QRCP、线性求解、SVD 等。`use_gpu=True` 在 GPU 实现前明确报错，不能像当前部分代码一样静默忽略。

## 2. 当前差距与处理顺序

现有 Python 已包含 `read/`、`gvector/`、`symmetry/`、`wavefunction/`、direct `epsilon.py`/`sigma.py` 和 HDF5 读写，但需先修正基线层：

1. `read/read_qe_gw.py` 的 `read_vxc` 尚未实现，k 点坐标/舍入、spinor 波函数轴序、占据数与 FFT 线性索引需逐项对齐。
2. `common/config.py` 的默认 `freq_dep`、`save_mem`、`no_symmetries_q_grid`、`exact_static_ch` 与 MATLAB 不一致，统一以 MATLAB 默认配置为准。
3. 当前 Python `sigma` 结果以 `eqp` 为主，尚未完整实现 MATLAB 的 `Vxc` 扣除、`eqp0`、`eqp1` 以及 `emf/ax/asx/ach/achx` 输出；先修复 direct，再引入 ISDF。
4. 现有 HDF5 restart 只恢复部分结果，不能替代 MATLAB 的兼容性检查与独立 sigma 恢复；统一为带 schema 和输入指纹的 checkpoint。
5. 当前测试主要覆盖 HDF5 round-trip 和少量 helper，缺少跨语言数值 golden 与端到端回归。

## 3. 目标结构与公共契约

在现有 `kssolv_spinor_gw/` 内增量整理为：

```text
common/       models、config、units、backend、progress、checkpoint 公共契约
read/         QE/KSSOLV/Vxc 输入
gvector/      G 矢量、cutoff、FFT 映射
symmetry/     full/irreducible BZ、subgroup、spinor symmetry
wavefunction/ genwf、FFT、normalization
coulG/        spherical 与 cell-box truncation
isdf/         options、sampling、space、solvers、polarizability、screened、cache
epsilon/      context、block preparation、direct/reduced operators、shared workflow
sigma/        context、block preparation、direct/reduced operators、QP、shared workflow
io/           HDF5 schema、checkpoint、qp.dat、golden 数据读写
cli.py        JSON 生产入口
```

关键决定：

- 保留 `Crystal`、`Syms` 等 dataclass；新增明确的 `EpsilonConfig`、`SigmaConfig`、`ISDFConfig` 与结果类型，边界处仍可提供 dict 兼容转换。
- 内部统一 0-based；仅在 MATLAB golden 导入/导出边界转换索引。能带轴固定为 `[band, k, spin]`，介电矩阵固定为 `[G, G', frequency]`，禁止各模块自行猜测轴序。
- 内部能量沿用 Ry，最终自能和 QP 输出为 eV；`ndiag_min/ndiag_max` 在 Python API 中为 0-based，JSON 提供显式 `index_base`，默认兼容 MATLAB 配置的 1-based 输入。
- `epsilon` 返回完整 screening 结果，可能是 `inv` 或 `isdf_screened_w`；`sigma` 根据 method 消费对应结果，不要求 reduced 路径额外构造完整 `eps.inv`。
- 配置校验集中完成，非法频率组合、band 范围、ISDF algorithm/solver/storage 组合立即报错。
- `common/backend.py` 只定义 `xp/fft/linalg/to_numpy` 等窄接口；本期实现 NumPy 后端，不把 CuPy 分支散落到数值公式中。

## 4. 分阶段实施

### 阶段 A：冻结参照与跨语言测试设施

- 建立功能—MATLAB 文件—Python 文件—测试用例对照表，参照 `src/GW/`、`test/ISDF/`、`test/cases/` 和 `test_mos2_gw_regression.m`。
- 在目标仓库增加 MATLAB golden 导出工具，将小型输入、固定插值点、中间数组和最终结果写成版本化 HDF5；记录 case id、配置、单位、轴序、索引基址和源代码 revision。
- 首批 golden：FFT/G 映射、对称性、genwf、direct matrix elements、chi0/epsilon、COHSEX/fullfreq sigma、Vxc/QP，以及固定点的 ISDF space/reduced W。
- 配置 `pytest`，将快速 unit、数据依赖 smoke、长时 validation 分开标记；保留现有 HDF5 测试并迁移到新 schema。
- 明确随机规则：Python 自身用固定 `numpy.random.Generator` 保证可复现；跨语言逐元素验证使用 MATLAB 导出的固定 `ind_mu`，不要求 MATLAB 与 NumPy 随机序列或 QRCP pivot 序列天然一致。

### 阶段 B：公共数值层与输入对齐

- 对齐 `read_qe_gw`：XML/HDF5 复数解码、spin/spinor 波函数、`mill`、k 点转换和六位舍入、密度、占据数；实现 `vxc.dat` 读取及 `read_vxc=False` 时的 `cal_vxc` 路径。
- 对齐 `gw_setup`、Fermi 能级、`nv`、G-vector/sort/gcutoff、FFT index、fullbz/irrbz/rqstar/subgrp/gmap、genwf 和 normalization。
- 补齐 Coulomb spherical/cell-box truncation，统一无截断与非法方法的处理。
- 用 scalar 和 spinor 小体系同时验证；发现 MATLAB 自身疑似缺陷时先记录差异并由参照测试确认，不在迁移中静默“修正”物理公式。

### 阶段 C：direct GW 完整对齐

- 将当前大型 `Epsilon`/`Sigma` 中的公共主循环与数值算子拆开，形成与 MATLAB `epsilon_ops`/`sigma_ops` 对应的 strategy 接口：`init/evaluate/accumulate/finalize` 和 `matrix_elements/contract`。
- direct epsilon 对齐 `save_mem` 两种累积方式、频率页、Coulomb 组装和矩阵求逆；求逆采用批量 `np.linalg.inv` 或等值 solve，不改变数学结果。
- direct sigma 对齐 q/k 映射、权重、SX/CH/X、exact static CH、contour deformation、频率存储与 `iw_lda`。
- 完成 `quasi_energy`、`get_eqp1`、symmetrize 和全部 QP 分量，统一输出 `sig/cor/emf/ax/asx/ach/achx/vxc/eqp0/eqp1`。
- 验收顺序：单个 matrix element → 单 q block → epsilon 全矩阵 → sigma 分量 → `eqp0/eqp1`，避免只用最终能量掩盖中间错误。

### 阶段 D：生产工作流、输出与恢复

- 新增 `run_json(config_path)` 和 `python -m kssolv_spinor_gw case.json`，字段与 MATLAB `gw_run_json.m` 对齐；所有相对路径相对于 JSON 文件，而不是当前工作目录。
- JSON 为主兼容入口，现有 YAML 示例保留为可选适配层，不形成第二套默认参数。
- 写出 BerkeleyGW 风格 `qp.dat`（含 `Eqp1`），结果 HDF5 保存配置、频率、G 矢量、screening 和 sigma/QP 数据。
- checkpoint 使用版本化 HDF5，保存 epsilon 请求配置、输入指纹、必要 system/symmetry 数据和结果；配置或 QE 输入不兼容时忽略缓存并重算，原子写入避免半成品。
- 支持仅加载 checkpoint 后运行 sigma，并覆盖 reduced `screened_w` 恢复，不强制生成完整逆矩阵。

### 阶段 E：ISDF 独立数值核

按 MATLAB `src/GW/ISDF/+isdf/` 的 public/private 边界迁移，Python 模块只暴露等价公共 API：

- `options.py`：rank/rank_ratio、product type、采样、求解器、stream 阈值、Cauchy、cache 参数的默认值与组合校验。
- `sampling.py`：QRCP（`scipy.linalg.qr(..., pivoting=True)`）、randomized QRCP、weighted/kmeans、unique fill、adaptive rank；支持 scalar 与 spinor component products。
- `space.py`：`build_space`、component products、product Gram、full/stream coefficients、`zeta_to_g`、real components、`matrix_elements`。
- `solvers.py`：direct stable solve、pinv 诊断、SVD whitening、条件数/rank 警告与 residual 检查。
- `polarizability.py`：direct 与 Cauchy reduced polarizability，保留 Frobenius 容差、迭代上限和系数分块策略。
- `screened.py`：`screened_w`、`screened_kernel`、`screened_contract`、full/paged solve，以及带字节预算和 hit/miss 统计的 screened-kernel cache。
- full 与 stream 只改变存储/分块，不改变结果；大 NN space、SVD 降秩和自适应 rank 均保留诊断元数据。

### 阶段 F：接入 epsilon/sigma 共享工作流

- method 统一解析为 `direct`、`matrix_elements`、`reduced_basis`；未提供或 `enable=False` 时必须走 direct，并忽略未启用的 ISDF 专用非法字段，与 MATLAB 行为一致。
- epsilon 接入 direct/matrix-elements/reduced accumulation、screened-W 输出及 irreducible/full-BZ 映射。
- sigma 接入 ISDF matrix elements、reduced screened contraction、VC/VN/NN product spaces、bare exchange space、全局 space 复用和 component/cache 生命周期。
- 覆盖 static 与 full-frequency reduced 路径；`screened_w` 默认模式验证不会创建 `eps.inv`。
- direct vs matrix-elements、direct vs full-inverse reduced、reduced `eps.inv` vs reduced screened-W 分别做等价性与物理误差验证。

### 阶段 G：性能、文档与发布验收

- 先保证正确，再优化：复用实空间波函数、减少临时数组、批量 BLAS/FFT、按元素数/字节数分块；full/stream 结果必须一致。
- 用 Si8/Si16 记录运行时间、峰值内存、ISDF rank、cache 命中率；大规模 MoS2 作为夜间/集群验证，不设未测量的任意加速比指标。
- 更新 `pyproject.toml`：加入 SciPy、pytest 开发依赖、CLI entry point 和包配置；修正当前 README 中与实际布局/启动命令不一致的内容。
- 文档提供 JSON 参数、Python API、checkpoint、CPU 后端限制、ISDF 调参和从 MATLAB 配置迁移示例。
- CPU 全功能验收后，再单独实现 CuPy backend 和 GPU 回归；不以 GPU 工作阻塞本期交付。

## 5. 验收矩阵

| 层级 | 用例 | 必须验证 |
|---|---|---|
| 合成单元测试 | scalar/spinor、单/多 component | FFT/G 映射、product cross terms、固定点插值、solver、screened kernel/cache |
| 小体系 smoke | Si2、Si8 direct/ISDF | JSON 全流程、exact static CH、qp.dat、checkpoint、full/stream 一致性 |
| 主要物理回归 | MoS2 2x2x2 | spinor 多 k COHSEX、direct fullfreq、ISDF reduced fullfreq、Eqp0/Eqp1 |
| 扩展验证 | AgBr、MoS2 4x4x4 | spinor ISDF/Cauchy、长时 fullfreq、内存与稳定性 |

- 确定性公共核和 direct 中间数组：与 MATLAB golden 使用 `allclose` 比较，常规目标 `rtol=1e-8, atol=1e-10`；求逆/病态求解按条件数单独设定并记录。
- MoS2 最终 `sig/eqp0/vxc/kpts` 沿用当前参照回归的严格标准：`abs_tol=1e-7`、`rel_tol=1e-10` 的 OR 判据。
- ISDF 固定 `ind_mu` 时比较插值基、矩阵元、reduced W 和 sigma；原生采样路径比较物理结果与误差趋势，不断言 pivot 索引完全一致。
- 所有缺数据或耗时用例必须显式标记并报告，不能用跳过测试替代核心 smoke 验收。

## 6. 完成定义

- 同一个 JSON 可在 MATLAB 与 Python 入口运行，Python 的路径解释、默认参数、输出字段和 `qp.dat` 与参照一致。
- direct、ISDF matrix-elements、ISDF reduced-basis 的 static/full-frequency、scalar/spinor、多 k 路径均有通过的核心回归。
- checkpoint 可跨进程恢复 direct 与 reduced sigma，且不会错误复用不兼容缓存。
- CPU 下无静默 GPU 参数、无未实现 `Vxc`/QP 分支、无仅靠 Python 自身结果生成的“伪基线”。
- 测试、示例、依赖和文档一致，后续 GPU 只需实现后端及性能验证，不重写数值工作流。
