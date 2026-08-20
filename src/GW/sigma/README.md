# Sigma Code Map

`sigma.m` owns the shared spin/band/k/q loop. `sigma_ops.m` chooses matrix
element construction and screened-interaction contraction handlers before
that loop.

## Handler Categories

Matrix element construction:

- `private/sigma_matrix_elements.m`

Full G-space contraction:

- `private/sigma_contract_full.m`

Reduced-basis contraction:

- `private/sigma_contract_reduced.m`

Shared output packing:

- `private/sigma_make_contribution.m`

Context, GPU, and mapping helpers:

- `private/sigma_kdata.m`
- `private/sigma_map_screened_w.m`
- `private/sigma_gpu_screened_w.m`
- `private/sigma_gather_if_gpu.m`
