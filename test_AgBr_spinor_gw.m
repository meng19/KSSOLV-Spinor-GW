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