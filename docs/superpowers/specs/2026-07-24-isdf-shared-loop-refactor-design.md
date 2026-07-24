# ISDF Shared-Loop Refactor Design

## 1. 背景

当前 `epsilon.m` 和 `sigma.m` 同时维护普通计算、ISDF matrix-elements 和 ISDF reduced-basis 三条路径。reduced-basis 通过在主函数顶部调用独立的大型工作流函数并提前 `return`，导致以下问题：

- k/q、自旋、波函数、FFT、对称映射和结果装配存在重复实现；
- GPU、频率模式和算法合法性检查散落在多个入口；
- `src/GW/ISDF` 下大量 `isdf_*` 函数平铺在全局 MATLAB 路径；
- 修改公共循环时必须同步多套工作流，容易产生行为差异；
- reduced-basis 的物理流程难以和普通参考路径逐块比较。

本次重构使用“共享主循环 + context struct + block struct + strategy function table”的折中方案。主函数继续完整展示物理循环，算法差异集中到少量策略接口，ISDF 数值算子迁移到 MATLAB `+isdf` package。

## 2. Git 检查点

重构前实现已经保存：

- commit: `e186aab` (`Checkpoint ISDF before shared-loop refactor`)
- tag: `isdf-before-shared-loop-refactor`

该 tag 是旧实现的唯一永久备份。源码中不长期保留 `epsilon_legacy.m`、`sigma_legacy.m` 或其他重复工作流。需要对照或恢复时使用 Git tag。

## 3. 目标

1. `epsilon.m` 只保留一套 q/spin/k 共享主循环。
2. `sigma.m` 只保留一套 spin/band/k/q 共享主循环。
3. 删除 reduced-basis 大函数调用后的提前 `return`。
4. direct、matrix-elements 和 reduced-basis 共用系统准备、波函数生成、FFT、对称映射、Coulomb 准备、星权重和结果存储。
5. 循环开始前只选择一次策略，循环内部不重复执行 `isfield`、`strcmpi` 或 `if use_isdf`。
6. 保留 reduced-basis 的低秩数据流；默认 `screened_w` 模式不得构造完整 `eps.inv`。
7. 使用 `+isdf` package 隔离 ISDF 数值算子并精简函数名。
8. 保持现有外部函数签名、配置字段和结果字段兼容。

## 4. 非目标

- 不在本次重构中增加 reduced-basis 全频计算。
- 不在本次重构中增加 reduced-basis GPU 支持。
- 不改变 ISDF 插值、Cauchy 求解、screened-W 或 COHSEX 的数学公式。
- 不为了形式统一而把 reduced-basis 强制展开成完整 G-space 矩阵。
- 不把所有数值叶子函数合并成单个巨型文件。

## 5. 外部兼容接口

外部调用保持不变：

```matlab
eps = epsilon(sys, options, syms, eps);
sig = sigma(eps, sig, sys, options, syms);
```

现有配置保持有效：

```text
eps.isdf.enable
eps.isdf.algorithm
eps.isdf.output
eps.isdf.reduced_solver
eps.isdf.sample_method

sig.isdf.enable
sig.isdf.algorithm
sig.isdf.reduced_solver
sig.isdf.sample_method
```

现有结果字段保持有效：

```text
eps.inv
eps.isdf_screened_w
eps.isdf_reduced_info
eps.isdf_reduced_rank

sig.sig
sig.cor
sig.eqp0
```

## 6. 核心数据模型

内部数据沿以下方向流动：

```text
context -> block -> contribution -> accumulator -> result
```

### 6.1 Context

`epsilon_context` 和 `sigma_context` 保存整次计算中不变或可复用的数据。context 在主循环前构造，之后按只读对象使用。

公共字段包括：

```text
sys, options, syms
method
use_gpu
nspin, nspinor, nbands
gvec, gr
normalization factors
FFT/grid metadata
irreducible/full-BZ mappings
Coulomb metadata
```

epsilon context 另外保存 `pol`、q 点 cutoff 和 k-star 映射。sigma context 另外保存 `fbz`、目标 k 点、epsilon inverse/full-BZ 映射和 screened-W/full-BZ 映射。

### 6.2 Block

block 是当前循环位置的临时输入，不保存跨 block 的可变状态。

epsilon block 至少包含：

```text
iq, ik, ispin
q, k, star maps/weight
wfnk, wfnkq
fft, idx
occupied bands, conduction bands
occupations and eigenvalues
```

sigma block 至少包含：

```text
ik, iq, in, ispin
k, q, q-star weight
wfnk, wfnkq
fft, idx
occupations
coulg, coulg_cutoff
eps_inv or mapped screened-W reference
exact-CH mapping metadata
```

### 6.3 Contribution

contribution 是策略对当前 block 的计算结果。

epsilon direct/matrix-elements contribution 保存完整 `chi0` block；epsilon reduced-basis contribution 保存低秩块：

```text
zeta_g
polar_coeff
rank
solver_info
```

sigma contribution 使用统一标量字段：

```text
asx
ax
ach
achx
```

### 6.4 Accumulator

accumulator 显式传入和返回，不依赖 persistent/global 隐藏状态。

- epsilon full-matrix accumulator 累加 `chi0`；
- epsilon reduced accumulator 保存 `zeta_g` 和 coefficient block 列表；
- sigma accumulator 累加 `asx/ax/ach/achx` 并应用 q-star 权重。

## 7. Strategy 接口

策略在进入循环前解析一次：

```matlab
ctx.method = resolve_method(input.isdf);
ops = epsilon_ops(ctx);
ops = sigma_ops(ctx);
```

`resolve_method` 只产生三种规范值：

```text
direct
matrix_elements
reduced_basis
```

### 7.1 Epsilon ops

```matlab
acc = ops.init(ctx, iq);
contribution = ops.evaluate(block);
acc = ops.accumulate(acc, contribution, block);
eps = ops.finalize(eps, acc, ctx, iq);
```

- `direct`: 使用普通矩阵元并形成完整 chi0 block；
- `matrix_elements`: 使用 ISDF 批量矩阵元并形成完整 chi0 block；
- `reduced_basis`: 构造 ISDF product space 和 reduced polarizability block，保持低秩累加。

`finalize` 负责 Coulomb 缩放、完整逆矩阵或 low-rank screened-W 的生成。`eps.isdf.output` 只影响 reduced-basis finalize，不改变共享循环。

### 7.2 Sigma ops

```matlab
matrix_elements = ops.matrix_elements(block);
contribution = ops.contract(block, matrix_elements);
```

- `direct`: 普通矩阵元 + 完整 `eps.inv` 收缩；
- `matrix_elements`: ISDF 批量矩阵元 + 完整 `eps.inv` 收缩；
- `reduced_basis`: ISDF product-space 矩阵元 + mapped screened-W 收缩；没有 screened-W 时允许按现有语义回退到完整 `eps.inv`。

裸交换、screened exchange、Coulomb hole 和 exact static CH 返回同一 contribution 结构，由主循环统一按 q-star 权重累加。

## 8. 共享主循环

### 8.1 Epsilon

```matlab
eps = epsilon_set_defaults(eps);
ctx = epsilon_context(sys, options, syms, eps);
ops = epsilon_ops(ctx);

for iq = 1:ctx.nq
    acc = ops.init(ctx, iq);
    for ispin = 1:ctx.nspin
        for ik = 1:ctx.nirred_k(iq)
            block = epsilon_prepare_block(ctx, iq, ik, ispin);
            contribution = ops.evaluate(block);
            acc = ops.accumulate(acc, contribution, block);
        end
    end
    eps = ops.finalize(eps, acc, ctx, iq);
end

eps.mtx = ctx.pol.mtx;
eps.nmtx = ctx.pol.nmtx;
eps.nfreq = ctx.pol.nfreq;
eps.freq = ctx.pol.freq;
```

q/k/spin、波函数、FFT、占据数和 k-star 映射始终在共享路径准备。direct 的能带对求和由其 block evaluator 完成；reduced-basis 对整个 occupied-conduction block 一次构造 product space。

### 8.2 Sigma

```matlab
sig = sigma_set_defaults(sig);
ctx = sigma_context(eps, sig, sys, options, syms);
ops = sigma_ops(ctx);

for ispin = 1:ctx.nspin
    for in = ctx.band_range
        for ik = 1:ctx.nk
            acc = sigma_init_accumulator();
            for iq = 1:ctx.nirred_q(ik)
                block = sigma_prepare_block(ctx, ik, iq, in, ispin);
                matrix_elements = ops.matrix_elements(block);
                contribution = ops.contract(block, matrix_elements);
                acc = sigma_accumulate(acc, contribution, block.weight);
            end
            sig = sigma_store(sig, acc, in, ik, ispin);
        end
    end
end

sig = sigma_finalize(sig, ctx);
```

主函数继续明确展示 `spin -> band -> k -> q` 物理顺序。策略调用替代散落的算法条件，但不会把主函数缩成只有一个 dispatcher 调用。

上面 `sigma_init_accumulator`、`sigma_accumulate`、`sigma_store` 和 `sigma_finalize` 表示主函数中现有的短小公共代码段，不新增同名工作流文件。它们负责零值初始化、四个标量贡献的加权相加、数组赋值以及最终 `quasi_energy`，不包含算法判断。实施时可以直接保留为主函数中的清晰代码；只有当某段需要独立单元测试时才提取为普通公共 helper。

## 9. 文件组织

```text
src/GW/epsilon/
  epsilon.m
  epsilon_context.m
  epsilon_prepare_block.m
  epsilon_ops.m

src/GW/sigma/
  sigma.m
  sigma_context.m
  sigma_prepare_block.m
  sigma_ops.m

src/GW/ISDF/+isdf/
  build_space.m
  polarizability.m
  screened_w.m
  screened_kernel.m
  screened_contract.m
  real_component.m
  matrix_elements.m

src/GW/ISDF/+isdf/private/
  product_gram.m
  component_products.m
  sample_points.m
  stable_solve.m
  zeta_to_g.m
```

工作流策略属于 epsilon/sigma，因此 `epsilon_ops.m` 和 `sigma_ops.m` 不放入 `+isdf`。`+isdf` 只包含真正的 ISDF 数值算子。

迁移后的主要名称为：

```text
isdf_build_space                 -> isdf.build_space
isdf_reduced_polarizability      -> isdf.polarizability
isdf_static_screened_interaction -> isdf.screened_w
isdf_screened_coulomb_kernel     -> isdf.screened_kernel
isdf_screened_coulomb_contract   -> isdf.screened_contract
isdf_wavefunction_real_component -> isdf.real_component
isdf_epsilon_batch               -> isdf.matrix_elements
isdf_sigma_batch                 -> isdf.matrix_elements
```

`isdf.matrix_elements` 是与 epsilon/sigma 无关的通用数值接口：

```matlab
gme = isdf.matrix_elements( ...
    left_components, right_components, idx_q, fftgrid, isdf_options);
```

epsilon/sigma strategy 负责从 block 选择左右波函数和能带，通用接口只接收 component product-space 数据。因此它不需要 `kind='epsilon'/'sigma'` 分支，也不依赖 GW 工作流结构。

数值叶子函数可以保留独立文件，但必须位于 package/private 边界内并只承担一个职责。不会用 operation-string dispatcher 把不相关算法塞入同一函数。

## 10. 错误处理

算法解析和模式约束在 context/ops 构造期间完成：

- 未知算法；
- reduced-basis 非静态模式；
- ISDF GPU 模式；
- 非法 epsilon output；
- sigma 缺少 screened-W 和 `eps.inv`；
- cutoff、G mapping、spinor component 或 band 维度不匹配。

进入循环后不再重复验证全局配置。block 相关错误应包含以下上下文：

```text
method, iq, ik, ispin, in
```

不使用宽泛 `try/catch` 隐藏底层异常；保留 MATLAB 原始调用栈和稳定的 `ISDF:*` error identifier。

## 11. 迁移顺序

### 阶段 1：行为刻画和公共数据结构

- 为 method resolution、context、block 和 ops selection 添加单元测试；
- 建立 `epsilon_context`、`sigma_context`；
- 建立 block preparation，但暂不删除旧路径；
- 确认公共准备结果与检查点实现一致。

### 阶段 2：共享 epsilon 循环

- 将 direct 和 matrix-elements 接入共享 epsilon loop；
- 将 reduced polarizability 作为 block contribution 接入；
- 将 full inverse/screened-W/both 移到 strategy finalize；
- 删除 `isdf_epsilon_reduced_basis(...); return;`；
- 完成 epsilon 单元和真实体系回归。

### 阶段 3：共享 sigma 循环

- 将 direct 和 matrix-elements 矩阵元构造接入 ops；
- 将 full inverse 和 mapped screened-W 收缩接入 ops；
- 统一 q-star 累加、exact static CH 和结果存储；
- 删除 `isdf_sigma_reduced_basis(...); return;`；
- 完成 sigma 单元和真实体系回归。

### 阶段 4：Package 迁移和清理

- 创建 `+isdf` 和 `+isdf/private`；
- 迁移并重命名数值算子；
- 更新内部调用和测试；
- 在迁移期保留最少量兼容 wrapper；
- 所有仓库内调用完成迁移后删除 wrapper 和旧大型工作流文件。

每个阶段独立通过测试后再进入下一阶段，不进行一次性全目录改写。

## 12. 测试设计

### 12.1 结构测试

- method resolution 对三种模式返回正确策略；
- context 必需字段和 mapping 尺寸正确；
- block preparation 对 scalar/spinor、单 K/多 K 正确；
- `epsilon.m` 和 `sigma.m` 不包含 reduced-basis 提前 return；
- 主循环中不重复解析 `isdf.enable/algorithm`；
- package private helper 不能作为公共 API 使用。

### 12.2 数值单元测试

- interpolation indices 和 sampling；
- scalar/spinor product space；
- spinor component cross terms；
- reduced polarizability direct/Cauchy；
- screened-W 与完整逆矩阵；
- epsilon/sigma matrix elements；
- irreducible/full-BZ G mapping；
- exact static CH。

### 12.3 工作流回归

- direct vs matrix-elements；
- direct vs reduced `full_inverse`；
- reduced `eps.inv` vs reduced screened-W sigma；
- `screened_w` 默认模式确认不创建 `eps.inv`。

### 12.4 真实体系

- AgBr：单 K、spinor、screened-W、exact static CH；
- MoS2 2x2x2：多 K、spinor、full-BZ q mapping、exact static CH。

所有数值结果以 tag `isdf-before-shared-loop-refactor` 和当前已通过测试为基线。重构不得放宽现有容差来掩盖差异。

## 13. 验收标准

1. `epsilon.m` 只有一套共享 q/spin/k 主循环。
2. `sigma.m` 只有一套共享 spin/band/k/q 主循环。
3. 两个主函数均无 ISDF 大工作流提前 return。
4. 算法只在循环前解析一次；循环内部通过 ops 调用差异算子。
5. direct、matrix-elements、reduced-basis 的现有输出字段兼容。
6. 默认 reduced epsilon 继续只保存 low-rank screened-W。
7. AgBr 和 MoS2 数值回归保持当前误差水平。
8. 所有 ISDF 单元测试、MATLAB Code Analyzer 和 `git diff --check` 通过。
9. `+isdf` package 外不再新增全局 `isdf_*` 内部辅助函数。
10. tag `isdf-before-shared-loop-refactor` 可以完整恢复重构前实现。
