# ISDF Cauchy COHSEX Design

## Goal

Add an optional cubic-scaling-oriented static COHSEX path based on ISDF and Cauchy-integral decoupling, following the Chapter 3 algorithm in Ma Huanhuan's thesis. Existing epsilon and sigma calculations remain the reference path and must keep their current defaults.

## Scope

The first implementation targets static COHSEX only:

- `eps.freq_dep = 0`
- `sig.freq_dep = 0`
- CPU execution
- Gamma or single-k validation first
- full-rank ISDF validation against the existing direct path

Full-frequency GW/CD, GPU execution, and full multi-k production support are outside the first implementation step. MoS2 multi-k validation can be added after the single-k path is verified.

## Architecture

The current ISDF path accelerates matrix-element construction but still feeds the existing high-scaling epsilon and sigma loops. The new path adds a separate low-rank static COHSEX route under `src/GW/ISDF/`, selected explicitly with `sig.isdf.algorithm = 'cauchy_cohsex'`.

The low-rank path builds separate ISDF spaces for the product pairs used in the thesis:

- `vc`: occupied by conduction states, used for polarizability and screened interaction.
- `vn`: occupied by target quasiparticle states, used for screened exchange.
- `nn`: summation states by target quasiparticle states, used for Coulomb-hole terms.

The Cauchy integral computes the reduced coefficient matrix `C * Omega^{-1} * C'` without explicitly looping over all occupied-unoccupied products in the final contraction. Screened interaction and self-energy are then formed in the reduced interpolation space.

## Components

New files:

- `src/GW/ISDF/isdf_build_space.m`: build interpolation indices, selected real-space values, and Fourier-space auxiliary basis for a pair of band sets.
- `src/GW/ISDF/isdf_cauchy_polarizability.m`: wrapper around `isdf_comega_cstar` for the `vc` coefficient matrix with direct and Cauchy modes.
- `src/GW/ISDF/isdf_static_screened_interaction.m`: construct the reduced static screened interaction from the `vc` basis and Coulomb kernel.
- `src/GW/ISDF/isdf_sigma_cohsex_cauchy.m`: compute static SEX and COH contributions in the reduced ISDF representation.

Modified files:

- `src/GW/sigma/sigma_set_defaults.m`: add `sig.isdf.algorithm`, `sig.isdf.cauchy_method`, `sig.isdf.cauchy_froErr`, and `sig.isdf.cauchy_MaxIter` defaults.
- `src/GW/sigma/sigma.m`: dispatch to the new static COHSEX path only when `sig.isdf.enable = true` and `sig.isdf.algorithm = 'cauchy_cohsex'`.

## Data Flow

1. `sigma.m` receives a normal static epsilon result and a sigma input with `sig.isdf.algorithm = 'cauchy_cohsex'`.
2. The new path validates that this is static COHSEX, CPU-only, and currently single-k compatible.
3. It builds `vc`, `vn`, and `nn` ISDF spaces from existing wavefunction objects and FFT helpers.
4. It computes the reduced polarizability coefficient matrix with `isdf_comega_cstar(..., method='cauchy')`.
5. It constructs reduced screened interaction matrices.
6. It computes SEX and COH contributions for the requested bands.
7. It returns a sigma-compatible output structure so existing tests and callers can compare `sig.sig` and `sig.eqp0`.

## Error Handling

Unsupported combinations fail explicitly instead of silently falling back:

- non-static frequency mode
- GPU execution
- missing band gap for Cauchy integration
- unsupported multi-k layout in the first implementation

The default calculation path remains unchanged when the new algorithm option is absent.

## Testing

Add tests under `test/ISDF/`:

- A reduced-operator unit test comparing direct and Cauchy `C * Omega^{-1} * C'`.
- An AgBr single-k static COHSEX validation comparing the new path with the existing direct static path at full ISDF rank.

Expected first-pass tolerance is relative error around `1e-6` to `1e-7`, tightened only after numerical behavior is confirmed.

## Acceptance Criteria

- Existing ISDF epsilon and sigma tests continue passing.
- The new AgBr Cauchy COHSEX validation passes against the existing direct static result.
- Default `epsilon` and `sigma` behavior is unchanged unless the new algorithm option is explicitly enabled.
- The implementation reduces the dominant polarizability contraction to the Cauchy-decoupled reduced-space form for the static COHSEX path.
