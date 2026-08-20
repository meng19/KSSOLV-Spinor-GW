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
