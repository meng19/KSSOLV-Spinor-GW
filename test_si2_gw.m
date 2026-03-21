% Cleaning up the workspace and command line window:
clc
clear all;
close all;
randn('state', 0);
rand('state', 0);
% Initializing the environment path for KSSOLV:
KSSOLV_startup;

read_vxc = 1;
[sys, options, syms] = read_qe_gw('.\example\qe_data\si2', read_vxc);
[sys, options] = gw_setup(sys, options);

tic
eps.nbnd = 10;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
eps.freq_dep = 2;
eps.freq_dep_method = 2;
eps.freq_cutoff = 200;
eps.delta_freq = 15;
eps.nfreq_imag = 15;
omega = 0;
eta = 0;
eps.cutoff = 5; %unit is Ry
% eps.coul_cut = 'spherical_truncation';
eps.coul_cut = 'cell_box_truncation';
eps.coul_cutoff = 5; %coulomb truncation radius in epsilon, only work with spherical_truncation
eps.use_gpu = 0;
eps.save_mem = 0;
eps = epsilon(sys, options, syms, eps);
toc

tic
sig.nbnd = 10;
sig.ndiag_min = 1;
sig.ndiag_max = 10;
% sig.coul_cut = 'spherical_truncation';
sig.freq_dep = 2;
sig.freq_dep_method = 2;
sig.freq_grid_shift = 2;
sig.max_freq_eval = 2;
sig.delta_freq_eval = 0.2;
sig.cd_int_method = 0;
sig.coul_cut = 'cell_box_truncation';
sig.coul_cutoff = 5; %coulomb truncation radius in sigma, only work with spherical_truncation
sig.no_symmetries_q_grid = 0;
sig.exact_static_ch = 0;
sig.use_gpu = 0;
sig = sigma(eps, sig, sys, options, syms);
toc