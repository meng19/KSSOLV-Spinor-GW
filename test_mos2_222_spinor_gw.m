% Cleaning up the workspace and command line window:
clc
clear all;
close all;
randn('state', 0);
rand('state', 0);
% Initializing the environment path for KSSOLV:
KSSOLV_startup;

read_vxc = 1;
[sys, options, syms] = read_qe_gw('.\example\qe_data\mos2_222_spinor', read_vxc);
[sys, options] = gw_setup(sys, options);

tic
eps.nbnd = 29;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
eps.freq_dep = 0;
% eps.freq_dep = 2;
% eps.freq_dep_method = 2;
% eps.freq_cutoff = 200;
% eps.delta_freq = 15;
% eps.nfreq_imag = 15;
omega = 0;
eta = 0;
eps.cutoff = 10;
eps.coul_cut = 'spherical_truncation';
eps.coul_cutoff = 10;
eps.use_gpu = 1;
eps.save_mem = 1;
eps = epsilon(sys, options, syms, eps);
% eps = epsilon_cohsex(sys, options, syms, eps);
toc

tic
sig.nbnd = 29;
sig.ndiag_min = 1;
sig.ndiag_max = 29;
sig.coul_cut = 'spherical_truncation';
sig.coul_cutoff = 10;
sig.no_symmetries_q_grid = 0;
sig.exact_static_ch = 1;
sig.use_gpu = 1;
sig = sigma(eps, sig, sys, options, syms);
toc