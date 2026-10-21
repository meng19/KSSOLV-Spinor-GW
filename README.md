# KSSOLV Spinor *GW*

[![License](https://img.shields.io/badge/License-[BSD_3_cluase]-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/Version-[2.0.1]-brightgreen.svg)](VERSION)

This is a spinor *GW* formalism within KSSOLV, a MATLAB toolbox for electronic structure calculations using Kohn-Sham density functional theory (DFT) which enables meticulous treatment of spinor in *GW* calculations.

## Table of Contents

- [KSSOLV Spinor *GW*](#kssolv-spinor-gw)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Installation](#installation)
  - [Directory Structure](#directory-structure)
  - [Example](#example)
  - [JSON configuration and parameter selection](#json-configuration-and-parameter-selection)
    - [Minimal workflow](#minimal-workflow)
    - [Common epsilon and sigma parameters](#common-epsilon-and-sigma-parameters)
    - [ISDF parameter guide](#isdf-parameter-guide)
    - [Recommended starting points](#recommended-starting-points)
  - [License](#license)
  - [How to cite](#how-to-cite)

## Features

- **New Function**: Supports spinor in *GW* calculations for molecules and periodic solids
- **Cross-platform**: Works on Windows, Linux, and macOS without compiling
- **Efficient**: Optimized for performance leveraging MATLAB's optimized linear algebra routines
- **High Performance**: Supports GPU parallelism and acceleration

## Installation

- **Prerequisites**: Ensure you have MATLAB installed (version R2019b or later recommended)
- **Verification**: Run a test script to verify the installation, e.g., `test_mos2_222_spinor_gw.m`
- **Run your file**: Similar to the example script, first edit your own script and save it in the folder root directory, then run it in MATLAB

## Directory Structure

The package is organized as follows:

```
kssolv-spinor-gw/
├── KSSOLV_startup.m     # Adds paths of the KSSOLV to Matlab
├── LICENSE              # BSD 3-Clause License text
├── README.md            # This file
├── src/                 # Main source code directory
│   ├── EigSolver/       # Diagonalization code
│   ├── GeomOpt/         # Structural optimization code
│   ├── GW/              # GW calculation code
        ├── common/      # General calculation functions, such as grid transformation, symmetry processing, etc.
        ├── epsilon/     # Calculate the dielectric matrix epsilon
        ├── kernel/      # Part to be improved, used for BSE calculation
        ├── read/        # Read wave functions, energy levels and other information for GW calculations
        ├── sigma/       # Calculate self-energy sigma
│   ├── SCF/             # KSSOLV ground state calculation
│   ├── Tools/           # Auxiliary functions such as visualization after ground state calculation
├── example/             # Example calculations files
│   ├── qe_data/         # The QE basis state calculation output results required in the example are stored here. Note that only HDF5 format results are compatible
│   └── ...
├── utils/               # Utility functions
├── ppdata/              # Pseudopotential files required for KSSOLV ground state calculations
└── external/            # External library files
```

## Example

- Simply run the `test_mos2_222_spinor_gw.m` file to calculate the quasi-partical energy of MoS<sub>2</sub> periodic solid or `test_AgBr_spinor_gw.m` for AgBr molecule

- Detailed description of `test_AgBr_spinor_gw.m`:

```MATLAB
% Cleaning up the workspace and command line window:
clc
clear all;
close all;
randn('state', 0);
rand('state', 0);
% Initializing the environment path for KSSOLV:
KSSOLV_startup;

% Whether to read the Vxc value of each band from Vxc.dat or to recalculate
read_vxc = 0;
% Read ground state wave function, energy level and other information from qe outputs
% Files needed: charge-density.hdf5, data-file-schema.xml, wfc*.hdf5(number of all k-points)
% Files optional: vxc.dat
[sys, options, syms] = read_qe_gw('.\example\qe_data\AgBr', read_vxc);
[sys, options] = gw_setup(sys, options);

% Epsilon calculation parameters
eps.nbnd = 30; % The number of energy bands in Epsilon calculation
eps.nv = options.nv; % Valence band number in Epsilon calculation
eps.nc = eps.nbnd - eps.nv; % Conduction band number in Epsilon calculation
eps.cutoff = 2; % Dielectric matrix cutoff in Epsilon calculations, in units of Ry
eps.coul_cutoff = 2; % Coulomb matrix cutoff in Epsilon calculations, in units of Bohr
eps.use_gpu = 0; % Whether to use GPU for Epsilon calculation
eps.save_mem = 0; % Whether to explicitly store the M matrix in the Epsilon calculation to speed up the summation of k-points and bands
eps = epsilon(sys, options, syms, eps); % Epsilon calculation main function

% Sigma calculation parameters
sig.nbnd = 30; % The number of energy bands in Sigma calculation
sig.ndiag_min = 1; % The lowest quasiparticle energy level number to be calculated in the Sigma calculation
sig.ndiag_max = 30; % The highest quasiparticle energy level number to be calculated in the Sigma calculation
sig.coul_cutoff = 2; % Coulomb matrix cutoff in Sigma calculations, in units of Bohr
sig.no_symmetries_q_grid = 0; % Whether k-point symmetry is considered in Sigma calculation
sig.exact_static_ch = 1; % Whether the static screened exchange is accurately calculated in the Sigma COHSEX calculation
sig.use_gpu = 0; % Whether to use GPU for Sigma calculation
sig = sigma(eps, sig, sys, options, syms); % Sigma calculation main function

% After the calculation is completed, the quasiparticle energy levels (ik, ib) of each k-point and band are stored in sig.eqp0.
```

## JSON configuration and parameter selection

`gw_run_json.m` is the recommended entry point for production runs.  It calls
`KSSOLV_startup` itself, so the submission script and the JSON input file may
live outside the repository.  Paths in the JSON file are interpreted relative
to the JSON file.  A complete, tested static Si example is
[`test/cases/si8_isdf.json`](test/cases/si8_isdf.json); the small smoke test is
[`test/cases/si8_isdf_smoke.json`](test/cases/si8_isdf_smoke.json).

### Minimal workflow

```matlab
addpath('C:/path/to/KSSOLV-Spinor-GW');
[eps, sig] = gw_run_json('C:/path/to/case.json');
```

The top-level JSON fields are:

| Field | Required | Meaning |
| --- | --- | --- |
| `qe_path` | yes | QE data directory.  The path is relative to this JSON file unless absolute. |
| `epsilon`, `sigma` | yes | Parameter objects passed to `epsilon` and `sigma`. |
| `read_vxc` | no, `true` | Read `vxc.dat` when it is available. |
| `rng_seed` | no | MATLAB random seed.  Fix it to make randomized ISDF point selection reproducible. |
| `qp_file` | no, `qp.dat` | Quasiparticle table written after sigma; relative to the JSON file. |
| `save_file` | no | MAT-file for `eps`, `sig`, system data, and the decoded configuration. |

### Common epsilon and sigma parameters

Set the physical convergence parameters before tuning ISDF.  In particular,
converge the number of bands and the dielectric cutoff against a direct (non-
ISDF) calculation first whenever feasible.

| Field | Where | Choice and effect |
| --- | --- | --- |
| `nbnd` | epsilon, sigma | Number of bands included.  Increase `epsilon.nbnd` for dielectric screening convergence; increase `sigma.nbnd` for the self-energy band sum. |
| `nv`, `nc` | epsilon | Explicit valence/conduction counts.  They default to the QE valence count and `nbnd-nv`; normally only `nbnd` is needed. |
| `ndiag_min`, `ndiag_max` | sigma | Inclusive band interval for quasiparticle corrections.  Reducing it saves sigma time but does not change the selected bands' physical band sum. |
| `cutoff` | epsilon | Dielectric-matrix kinetic-energy cutoff (Ry).  A main physical convergence parameter; larger values raise memory and time. |
| `coul_cut`, `coul_cutoff` | epsilon, sigma | Coulomb treatment and its cutoff.  Use matching values in epsilon and sigma.  For example, the Si cases use `"spherical_truncation"` and `5`. |
| `freq_dep` | epsilon, sigma | `0` selects static COHSEX.  `2` enables the frequency-dependent workflow; retain the corresponding frequency-grid inputs used by that workflow. |
| `use_gpu` | epsilon, sigma | Move supported array work to a MATLAB GPU.  Keep `false` if GPU memory is insufficient. |
| `precompute_wav` | epsilon, sigma | Cache wavefunctions for less repeated I/O/FFT work at the cost of memory. |
| `save_mem` | epsilon | Store less intermediate wavefunction/matrix-element data; trades memory for extra computation. |
| `no_symmetries_q_grid` | sigma | `false` uses q-grid symmetries (recommended when valid); `true` disables symmetry reduction for debugging or a deliberately full grid. |
| `exact_static_ch` | sigma | Use the direct static Coulomb-hole expression instead of the default reduced/static path.  Use as a validation reference; it is generally more expensive. |

### ISDF parameter guide

ISDF is activated independently in `epsilon.isdf` and `sigma.isdf`.  The
static reduced-basis configuration below is the normal high-performance
choice.  `matrix_elements` retains ISDF matrix-element construction but uses
the conventional full-space downstream workflow; it is useful for validation
and frequency-dependent runs, not the usual memory-saving option.

```json
"isdf": {
  "enable": true,
  "algorithm": "reduced_basis",
  "sample_method": "qrcp_randomized",
  "sample_precision": "double",
  "interpolation_solver": "svd_whiten",
  "svd_cutoff": 1e-12,
  "svd_ratio": 0.75,
  "seed": 0
}
```

| Field | Valid values / default | How to choose |
| --- | --- | --- |
| `enable` | `false` | Set `true` only after a direct reference is available for the system. |
| `algorithm` | `reduced_basis` or `matrix_elements` | Use `reduced_basis` for static ISDF screening/self-energy.  Use `matrix_elements` to validate matrix elements or for a workflow that still needs the full dielectric representation. |
| `rank`, `rank_ratio` | positive; ratio default `1` | ISDF rank is `ceil(sqrt(nleft*nright)*rank_ratio)`, capped by the grid/product-space size.  Use an explicit `rank` only for a convergence test; otherwise use a ratio.  Larger ranks improve accuracy but increase the cubic reduced-space operations. |
| `rank_ratio_vc` | epsilon override | Rank ratio for the valence--conduction (VC) space used to construct screening.  It takes precedence over the generic `rank_ratio`. |
| `rank_ratio_nn` | sigma override | Rank ratio for the target-band--all-band (NN) space used in screened exchange/COH sigma contractions. |
| `rank_ratio_vn` | sigma override | Rank ratio for the target-band--occupied-band (VN) space used by bare exchange when selected. |
| `sample_method` | `qrcp` (default), `qrcp_randomized`, `kmeans` | `qrcp_randomized` is the recommended speed/quality choice for large scalar product spaces.  `qrcp` materializes the product matrix and is a more expensive deterministic reference.  `kmeans` is a cheaper alternative that should be convergence-tested. |
| `seed` | `0` | Fix it for repeatability.  Different seeds can select different interpolation points and give slightly different compressed results. |
| `sample_precision` | `double` (default), `single` | `single` applies only to temporary randomized sampling sketches and can speed the QRCP stage/reduce its memory.  Keep `double` for production accuracy checks; validate before using `single`. |
| `projection_block_elements` | `Inf` | The default uses one complete projected-product block for maximum throughput.  Set a finite limit only for memory-limited calculations. |
| `rank`, `rank_ratio`, `warn_rank_selection` | recommended rank, `1`, `true` | `rank` defaults to `ceil(sqrt(nleft*nright)*rank_ratio)` after clamping to the product-space size.  The warning reports deliberate over/under-ranking. |
| `weight`, `power` | `add`, `1` | K-means point weight.  For an implicit multi-component weight, the implementation uses `prod`; explicit `add`, `prod`, or `power` overrides it. |
| `init`, `kmeans_max_iter`, `kmeans_distance_block_elements` | `random`, `100`, `Inf` | K-means initializer, iteration cap, and maximum temporary distance-matrix size.  Use a finite element limit only to bound memory. |
| `random_oversampling` | `1.2` | Oversampling factor for randomized QRCP sketches. |
| `product_mu_storage`, `product_mu_max_elements` | `full`, `1e8` | Band-pair coefficient representation.  Select `auto` or `stream` explicitly for memory-limited calculations. |
| `product_gram_storage`, `product_gram_max_elements`, `product_gram_block_elements` | `full`, `1e8`, `Inf` | Gram/interpolation representation.  `stream` uses a `2e6`-element block unless a finite block limit is supplied. |
| `rcond_tol`, `warn_ill_conditioned`, `ill_conditioned_solver` | `1e-12`, `false`, `direct` | Conditioning policy for interpolation solves.  Set solver to `pinv` only for a deliberate pseudoinverse fallback. |
| `swap_left_right` | `false` | Validation control that selects points from the conjugate right-left space while retaining the original physical representation. |
| `rank_vc`, `rank_vn`, `rank_nn`; `rank_ratio_vc`, `rank_ratio_vn`, `rank_ratio_nn`; `swap_left_right_vc`, `swap_left_right_vn`, `swap_left_right_nn` | unset | Optional type-specific overrides.  When unset, the corresponding base `rank`, `rank_ratio`, or `swap_left_right` value is used. |
| `points`, `fftgrid` | `[]`, build input | Optional K-means coordinates and the real-space FFT grid.  `fftgrid` is populated from `isdf.build_space` when omitted. |
| `progress`, `print_rank` | `[]`, `true` | Optional progress callback and rank-report control. |
| `validate_hf_exchange_max_pairs`, `validate_hf_exchange_label` | `Inf`, `HF exchange` | Limit and label for the optional direct-vs-ISDF HF exchange diagnostic.  Epsilon/sigma own the separate boolean that enables this diagnostic. |

| `interpolation_solver` | `direct` (default), `svd_whiten` | `direct` uses the original interpolation coordinates.  `svd_whiten` applies the product-Gram SVD normalization and can truncate ill-conditioned modes; it is recommended for source-compatible static ISDF runs. |
| `svd_cutoff` | `0` | Relative eigenvalue threshold for `svd_whiten`; modes below `max_eigenvalue*svd_cutoff` are removed.  Use `0` for a pure coordinate transform, or start with `1e-12` when conditioning is a concern.  Any truncation must be convergence-tested. |
| `svd_ratio` | `0.5` | Split of the Gram inverse between interpolation basis and coefficients.  It has no physical effect without truncation.  `0.75` reproduces the convention used by the reference ISDF case. |
| `adaptive_rank_enable` | `false` | Increase rank until an interpolation residual target is reached.  Requires QRCP sampling and is incompatible with `svd_whiten`; use it only for rank exploration. |
| `adaptive_rank_tol`, `adaptive_rank_step`, `adaptive_rank_max`, `adaptive_validation_rank` | `1e-6`, `ceil(sqrt(nleft*nright))`, maximum allowed rank, `max(16, adaptive_rank_step)` | Controls adaptive-rank tolerance, increment, maximum rank, and randomized validation sketch.  Leave unset unless `adaptive_rank_enable=true`. |
| `reduced_solver` | `cauchy` (default), `direct` | Solver for the reduced polarizability.  `cauchy` is the efficient static approximation/iteration and falls back to direct pages when needed; `direct` is the reference sum and is useful for validating Cauchy settings. |
| `cauchy_froErr`, `cauchy_MaxIter` | `1e-8`, `12` | Requested Cauchy relative Frobenius error and iteration cap.  Tighten/increase only if its fallback summary or a direct comparison shows a need. |
| `cauchy_coefficient_storage` | `full` (default), `auto`, `stream` | Controls whether Cauchy pair coefficients are cached.  The default keeps all P/Q coefficients and uses one block for maximum throughput.  Use `auto` or `stream` for memory-limited large systems. |
| `cauchy_coefficient_cache_bytes`, `cauchy_coefficient_block_bytes` | `4 GiB`, `Inf` for `full`; `512 MiB` block for `auto`/`stream` | Byte limits for the full P/Q cache and a streamed P/Q block.  Set finite limits only when memory must be bounded. |

All ISDF storage controls use the same mode names: `full` materializes the complete object, `stream` creates bounded temporary blocks, and `auto` selects between them from its size limit.

`isdf.workflow_defaults` is the single source for epsilon/sigma shared workflow defaults.  Product-space and numerical defaults are centralized in `isdf/private/set_defaults`, and are applied when a product space or reduced polarizability is actually built.  `max_rank`, `recommended_rank`, `unclamped_recommended_rank`, `rank_source`, and `weight_was_set` are internal derived fields populated by `set_defaults`; `product_type` is internal context populated by `isdf.options_for_type`.  They are not JSON configuration parameters.  Likewise, the internal polarizability fields `method`, `freq`, `froErr`, and `MaxIter` are populated centrally by `set_defaults` after the epsilon/sigma configuration fields `reduced_solver`, `cauchy_froErr`, and `cauchy_MaxIter` are mapped into the solver input.

`epsilon.isdf.output` applies only to `algorithm="reduced_basis"`:

| Value | Result and use |
| --- | --- |
| `"screened_w"` (recommended) | Keep only the reduced screened interaction needed by ISDF sigma; lowest memory. |
| `"full_inverse"` | Build only the full reciprocal-space dielectric inverse; use for compatibility/diagnostics. |
| `"both"` | Keep both representations; highest memory. |

It is an output-storage choice, not a switch between SMW and direct matrix
inversion.  The reduced screened-W construction is used in the reduced-basis
path in all three modes.

Sigma has further product-space reuse controls:

| Field | Default | Meaning and recommendation |
| --- | --- | --- |
| `exchange_space` | `"nn"` | Product space for bare exchange: `"nn"` reuses the NN space; `"vn"` builds/uses a smaller occupied-band VN space.  `"vn"` can be faster for few occupied bands, but needs its own rank convergence. |
| `reuse_nn_for_vn` | `false` | If `true`, force exchange to reuse NN even when `exchange_space="vn"`; no VN space is built.  This prioritizes reuse over the smaller VN representation. |
| `global_nn_space` | `false` | If `true`, build one NN space for every requested diagonal band at a given k/q/spin, then reuse it.  It usually reduces repeated setup when several `ndiag` bands are requested, but may use more memory. |
| `global_vn_space` | `false` | Analogous reuse for the VN exchange space.  It only matters when a separate VN space is actually used. |
| `reuse_eps_real_wfn` | `false` | Reuse epsilon's cached real-space wavefunctions in sigma.  Set it together with `epsilon.isdf.cache_real_wfn=true`; it saves FFT construction at the cost of retaining the cache. |
| `cache_real_wfn` | `false` (epsilon) | Retain epsilon real-space wavefunction components so sigma can reuse them.  Helpful when epsilon and sigma use compatible FFT grids; consumes memory. |
| `screened_kernel_cache_bytes` | `Inf` | Byte budget for reusing a projected reduced screened kernel across the diagonal bands that share one target product space.  Reuse pays off with `global_nn_space=true`; the exact static CH matrix depends only on q and is reused regardless.  `Inf` (default) keeps every reusable kernel without a memory cap; entries are never evicted, so set a finite byte cap — or `0` to disable reuse — when memory is tight.

The progress text reflects these choices: `VC`, `NN`, and `VN` name the
product spaces; `NN built: 32 diag / 319 sum` means an NN space was
constructed for 32 requested diagonal bands and 319 summation bands; `NN
reused` means the previously built global NN space was reused.  It is
unrelated to `reuse_nn_for_vn`.

`screened_kernel_cache_bytes` covers the third setup cost of a reduced-basis
sigma run: projecting the reduced screened interaction onto a target product
space.  With `global_nn_space=true` one projection serves every requested
diagonal band, so the run builds it once per k/q instead of once per band; the
exact static CH matrix is q-only and is reused in either mode.  The
`Sigma screened kernel` row of the timing report counts those projections and
shows the miss cost in `CPU max`, and the run ends with an
`ISDF screened-kernel cache: ...` summary line.  Per-band NN spaces
(`global_nn_space=false`) cannot share a projection, so the budget only helps
when a global space or the exact CH term is used.

### Recommended starting points

1. **Physical baseline:** run direct (`isdf.enable=false`) with converged
   `epsilon.nbnd`, `sigma.nbnd`, `epsilon.cutoff`, and Coulomb settings.
2. **Static production ISDF:** start from `test/cases/si8_isdf.json`:
   `reduced_basis`, `qrcp_randomized`, `double`, `svd_whiten`,
   `svd_cutoff=1e-12`, `svd_ratio=0.75`, `reduced_solver="cauchy"`, and fixed
   `seed=0`.  Raise VC, NN, and (if used) VN rank ratios separately until
   `SX-X`, `CH`, and `Eqp0` agree with the direct baseline at the desired
   tolerance.
3. **Fast exploratory calculation:** lower the three rank ratios and consider
   `sample_precision="single"`; compare the resulting `qp.dat` with the
   double-precision production case before trusting it.
4. **Many QP bands:** enable `global_nn_space`; if using a distinct VN exchange
   space, also enable `global_vn_space`.  Monitor memory, since the global
   spaces retain all requested diagonal bands together.  The global space also
   lets `screened_kernel_cache_bytes` (default `Inf`, unlimited) replace
   `ndiag - 1` repeated screened projections per k/q by one; set a finite
   budget there if the retained kernels grow too large.
5. **Diagnosing a discrepancy:** set `reduced_solver="direct"`, keep
   `sample_precision="double"`, use a fixed seed, and compare one product
   space/rank at a time against the direct calculation.  Avoid simultaneously
   changing rank, sampling method, and SVD cutoff.

## License

This software is licensed under the BSD 3-Clause License, one of the more permissive free software licenses. This license allows you to use, modify, and distribute the software either in source code or binary form. For the specific terms and conditions, refer to the full license text available online at [https://opensource.org/license/BSD-3-Clause](https://opensource.org/license/BSD-3-Clause)

## How to cite
