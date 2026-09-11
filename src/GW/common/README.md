# GW Common Utilities

This directory contains public GW helper functions shared by epsilon, sigma,
the input readers, and Coulomb kernels.  The subdirectories are grouped by
numerical role; `KSSOLV_startup` adds `src/` recursively, so these moves do not
change MATLAB function names or call sites.

| Directory | Contents |
| --- | --- |
| `reciprocal/` | G-vector construction, reciprocal/grid mappings, cutoff selection, FFT-grid bounds, and reciprocal-vector sorting. |
| `symmetry/` | k/q-point matching, irreducible/full Brillouin-zone construction, and symmetry-group helpers. |
| `electronic/` | Occupations, Fermi level, and exchange-correlation potential helpers. |
| `wavefunction/` | Wavefunction generation and normalization checks. |
| `runtime/` | GW setup, workflow selection, progress/timing/reporting, GPU gathering, and FFT-size checks. |

## Deliberately unchanged areas

- `epsilon/private/` and `sigma/private/` remain flat because MATLAB resolves
  private functions only in their parent directory.  Adding another directory
  level would either break that visibility rule or turn internal helpers into
  globally visible functions.
- `ISDF/+isdf/private/` is likewise left intact: it is the private scope of a
  MATLAB package.  Its file names and the package README already group the
  code into sampling, product-space construction, reduced polarizability, and
  screened-interaction helpers.
- `read/`, `coulG/`, and `kernel/` are each small, cohesive functional units;
  additional nesting would not make navigation easier.
