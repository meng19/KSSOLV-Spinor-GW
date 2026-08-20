# Epsilon Code Map

`epsilon.m` owns the shared q/spin/k loop. `epsilon_ops.m` chooses the
workflow-specific handlers once before that loop.

## Handler Categories

Full G-space path:

- `private/epsilon_full_init.m`
- `private/epsilon_full_evaluate.m`
- `private/epsilon_full_accumulate.m`
- `private/epsilon_direct_stream_accumulate.m`
- `private/epsilon_full_finalize.m`

Full matrix assembly helpers:

- `private/epsilon_add_state_batch.m`
- `private/epsilon_mapped_gme.m`
- `private/epsilon_repeat_eden.m`
- `private/epsilon_invert_pages.m`

Context and loop helpers:

- `private/epsilon_qdata.m`
- `private/epsilon_prepared_data.m`

Reduced-basis path:

- `private/epsilon_reduced_init.m`
- `private/epsilon_reduced_evaluate.m`
- `private/epsilon_reduced_accumulate.m`
- `private/epsilon_reduced_finalize.m`

Reduced-basis assembly helpers:

- `private/epsilon_reduced_accumulate_mapped_chi.m`
- `private/epsilon_reduced_mapped_zeta_chi.m`
- `private/epsilon_page_blkdiag.m`
- `private/epsilon_repeat_page_blkdiag.m`
- `private/epsilon_gather_screened_w.m`
