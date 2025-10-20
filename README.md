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

## License

This software is licensed under the BSD 3-Clause License, one of the more permissive free software licenses. This license allows you to use, modify, and distribute the software either in source code or binary form. For the specific terms and conditions, refer to the full license text available online at [https://opensource.org/license/BSD-3-Clause](https://opensource.org/license/BSD-3-Clause)

## How to cite
