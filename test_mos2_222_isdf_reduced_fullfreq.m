% MoS2 2x2x2 multi-k spinor full-frequency GW calculation with ISDF only.
clc;
clear;
close all;
rng(0, 'twister');

script_dir = fileparts(mfilename('fullpath'));
repo_root = script_dir;
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

[sys, options, syms] = read_qe_gw( ...
    '.\example\qe_data\mos2_222_spinor', 1);
[sys, options] = gw_setup(sys, options);
assert(sys.nkpts > 1, 'This script requires a multi-k data set.');
assert(sys.nspinor > 1, 'This script requires spinor wavefunctions.');

tic
eps.nbnd = 29;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
eps.freq_dep = 2;
eps.freq_dep_method = 2;
eps.freq_cutoff = 200;
eps.delta_freq = 15;
eps.nfreq_imag = 15;
eps.cutoff = 10;
eps.coul_cut = 'spherical_truncation';
eps.coul_cutoff = 10;
eps.use_gpu = 0;
eps.save_mem = 0;
eps.precompute_wav = 0;
eps.isdf.enable = true;
eps.isdf.algorithm = 'reduced_basis';
eps.isdf.output = 'screened_w';
eps.isdf.sample_method = 'kmeans';
eps.isdf.rank = 28;
eps.isdf.reduced_solver = 'direct';
eps.isdf.seed = 0;
eps = epsilon(sys, options, syms, eps);
fprintf('Finished ISDF reduced-basis epsilon screened W.\n');
toc

tic
sig.nbnd = 29;
sig.ndiag_min = sig.nbnd;
sig.ndiag_max = sig.nbnd;
sig.freq_dep = 2;
sig.freq_dep_method = 2;
sig.freq_grid_shift = 2;
sig.max_freq_eval = 2;
sig.delta_freq_eval = 0.2;
sig.cd_int_method = 0;
sig.coul_cut = 'spherical_truncation';
sig.coul_cutoff = 10;
sig.no_symmetries_q_grid = 1;
sig.exact_static_ch = 1;
sig.use_gpu = 0;
sig.precompute_wav = 0;
sig.isdf.enable = true;
sig.isdf.algorithm = 'reduced_basis';
sig.isdf.sample_method = 'kmeans';
sig.isdf.rank = 28;
sig.isdf.seed = 0;
sig = sigma(eps, sig, sys, options, syms);
fprintf('Finished ISDF reduced-basis sigma full-frequency calculation.\n');
toc
