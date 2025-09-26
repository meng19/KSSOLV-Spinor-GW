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
omega = 0;
eta = 0;
eps.cutoff = 10; %unit is Ry
eps.coul_cutoff = 10; %coulomb truncation radius in epsilon
eps.use_gpu = 1;
eps.save_mem = 1;
eps = epsilon(sys, options, syms, eps);
toc

tic
sig.nbnd = 29;
sig.ndiag_min = 1;
sig.ndiag_max = 29;
sig.coul_cutoff = 10; %coulomb truncation radius in sigma
sig.no_symmetries_q_grid = 0;
sig.exact_static_ch = 1;
sig.use_gpu = 1;
sig = sigma(eps, sig, sys, options, syms);
toc