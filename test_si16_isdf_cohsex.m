% Cleaning up the workspace and command line window:
clc
clear all;
close all;
randn('state', 0);
rand('state', 0);
isdf = true;

% 启动 Profiler
% profile on -memory

% Initializing the environment path for KSSOLV:
KSSOLV_startup;

read_vxc = 1;
[sys, options, syms] = read_qe_gw('.\example\qe_data\si16', read_vxc);
[sys, options] = gw_setup(sys, options);

eps.nbnd = 127;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
eps.freq_dep = 0;
eps.cutoff = 5; % Ry
eps.coul_cut = 'spherical_truncation';
eps.coul_cutoff = 5; %coulomb truncation radius in epsilon
eps.save_mem = 0;
eps.use_gpu = 0;
eps.precompute_wav = 0;
eps.isdf.enable = isdf;
eps.isdf.algorithm = 'reduced_basis';
eps.isdf.output = 'screened_w';
eps.isdf.sample_method = 'qrcp_randomized';
%eps.isdf.rank_ratio = 16;
eps.isdf.rank = 64;
eps.isdf.reduced_solver = 'cauchy';
eps.isdf.seed = 0;

tic
eps = epsilon(sys, options, syms, eps);
toc

sig.nbnd = 127;
sig.ndiag_min = 120;
sig.ndiag_max = 127;
sig.freq_dep = 0;
sig.coul_cut = 'spherical_truncation';
sig.coul_cutoff = 5; %coulomb truncation radius in sigma
sig.no_symmetries_q_grid = 0;
sig.exact_static_ch = 1;
sig.use_gpu = 0;
sig.precompute_wav = 0;
sig.isdf.enable = isdf;
sig.isdf.algorithm = 'reduced_basis';
sig.isdf.sample_method = 'qrcp';
%sig.isdf.rank_ratio = 16;
sig.isdf.rank = 12;
sig.isdf.seed = 0;

tic
sig = sigma(eps, sig, sys, options, syms);
toc

% fid = fopen('eqp0.txt', 'w');
% fprintf(fid, [repmat('%.15g\t', 1, size(sig.eqp0, 2)-1), '%.15g\n'], sig.eqp0.');
% fclose(fid);
% 
% profile off
% p = profile('info');
% profsave(p, 'profile');

