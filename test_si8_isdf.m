% Si8 Gamma COHSEX/ISDF case imported from GW_submodule/examples/cases.
% The source namelist uses vc/vn/nn slots.  The current implementation
% maps vc to epsilon and nn to sigma; it has no separate vn slot yet.

clc;
clear;
close all;
rng(0, 'twister');

KSSOLV_startup;

read_vxc = 1;
[sys, options, syms] = read_qe_gw('.\\example\\qe_data\\si8', read_vxc);
[sys, options] = gw_setup(sys, options);

% External case: number_bands_in_summation = 319, coulomb_cutoff = 15,
% spherical truncation radius = 5, static Gamma COHSEX.
eps = struct();
eps.nbnd = 319;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
eps.freq_dep = 0;
eps.cutoff = 15;
eps.coul_cut = 'spherical_truncation';
eps.coul_cutoff = 5;
eps.save_mem = false;
eps.use_gpu = false;
eps.precompute_wav = false;
eps.isdf.enable = true;
eps.isdf.algorithm = 'reduced_basis';
eps.isdf.output = 'screened_w';
% The reference repository's exxmethod='qrcp' is randomized QRCP.
eps.isdf.sample_method = 'qrcp_randomized';
eps.isdf.rank_ratio_vc = 12;
eps.isdf.random_oversampling = 1.2;
% In the source case QRCP is a final-point selector; its adaptive fields
% apply only to the source repository's pseudo backend.
eps.isdf.adaptive_rank_enable = false;
eps.isdf.seed = 0;
eps.isdf.reduced_solver = 'cauchy';

fprintf('Running imported Si8 epsilon ISDF case...\n');
tic;
eps = epsilon(sys, options, syms, eps);
toc;

sig = struct();
sig.nbnd = 319;
sig.ndiag_min = 1;
sig.ndiag_max = 32;
sig.freq_dep = 0;
sig.coul_cut = 'spherical_truncation';
sig.coul_cutoff = 5;
sig.no_symmetries_q_grid = false;
sig.exact_static_ch = false;
sig.use_gpu = false;
sig.precompute_wav = false;
sig.isdf.enable = true;
sig.isdf.algorithm = 'reduced_basis';
sig.isdf.sample_method = 'qrcp_randomized';
sig.isdf.random_oversampling = 1.2;

sig.isdf.adaptive_rank_enable = false;
% Build NN once for all diagonal bands and reuse it, matching the source
% repository's type-space lifetime instead of rebuilding once per band.
sig.isdf.rank_ratio_nn = 16;
sig.isdf.rank_ratio_vn = 24;
sig.isdf.exchange_space = 'vn';
sig.isdf.reuse_nn_for_vn = true;
sig.isdf.global_nn_space = true;
sig.isdf.global_vn_space = true;
sig.isdf.seed = 0;

fprintf('Running imported Si8 sigma ISDF case...\n');
tic;
sig = sigma(eps, sig, sys, options, syms);
toc;

fprintf('Imported Si8 ISDF case completed.\n');
