# ISDF Package Code Map

This package contains the reusable ISDF numerical kernels used by GW
epsilon and sigma workflows.

## Public Package API

- `build_space.m`: build an ISDF product space from real-space wavefunction
  components, selected interpolation points, and reciprocal-space indices.
- `matrix_elements.m`: build Fourier matrix elements using an ISDF product
  space.
- `polarizability.m`: build reduced polarizability coefficients. Solver
  details live in `private/comega_*.m`.
- `real_component.m`: extract one real-space wavefunction component.
- `screened_w.m`: construct the reduced screened interaction object.
- `screened_kernel.m`: project reduced screened interaction into a target
  product space.
- `screened_contract.m`: contract a reduced screened kernel with product
  coefficients.

## Private Helper Categories

Sampling:

- `sample_points.m`: dispatch by `sample_method`.
- `qrcp_sample.m`: QRCP point selection.
- `scalar_randomized_sample.m`: randomized QRCP for scalar products.
- `kmeans_sample.m`, `kmeanspp_init.m`, `grid_points.m`,
  `distance_to_centers.m`: weighted K-means selection.
- `weighted_sample.m`, `weighted_one.m`, `unique_fill.m`: sampling repair
  and weighted draws.

Product-space construction:

- `component_products.m`: explicit component product matrix construction.
- `component_weight.m`: O(ngrid) component-product weights for K-means.
- `pair_products.m`: scalar pair-product materialization.
- `product_gram.m`: Gram matrices for the interpolation solve.
- `set_defaults.m`: product-space option defaults.
- `zeta_to_g.m`: transform real-space interpolation basis to reciprocal
  space.
- `stable_solve.m`: guarded right solve for interpolation systems.

Reduced polarizability:

- `comega_cstar.m`: solver dispatch.
- `comega_direct.m`: explicit occupied-unoccupied summation.
- `comega_cauchy.m`: static Cauchy solver.

Common utilities:

- `gather_if_gpu.m`: gather `gpuArray` values only when needed.
- `page_project_kernel.m`: page-wise screened-kernel projection.
- `page_solve_screened_w.m`: page-wise reduced screened-W linear solves.
- `randn_like.m`: random normal values matching reference precision/device.

## Review Notes

The package API is intentionally small. New external calls should prefer the
public functions above; helper files under `private/` are package-internal
implementation details. When reviewing changes, start from the public entry
point, then follow only the helper category used by the changed option or
workflow.
