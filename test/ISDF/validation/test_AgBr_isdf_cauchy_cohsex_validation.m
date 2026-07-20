% Validate spinor ISDF Cauchy COHSEX sigma path against existing direct path.
clc;
clear all;
close all;
randn('state', 0);
rand('state', 0);

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

read_vxc = 0;
[sys, options, syms] = read_qe_gw('.\example\qe_data\AgBr', read_vxc);
[sys, options] = gw_setup(sys, options);

eps_input.nbnd = options.nv + 2;
eps_input.nv = options.nv;
eps_input.nc = eps_input.nbnd - eps_input.nv;
eps_input.cutoff = 2;
eps_input.coul_cut = 'cell_box_truncation';
eps_input.use_gpu = 0;
eps_input.save_mem = 1;
eps_input.precompute_wav = 0;
eps_input.freq_dep = 0;

eps_result = epsilon(sys, options, syms, eps_input);

sig_base.nbnd = eps_input.nbnd;
sig_base.ndiag_min = options.nv;
sig_base.ndiag_max = options.nv;
sig_base.coul_cut = 'cell_box_truncation';
sig_base.no_symmetries_q_grid = 0;
sig_base.exact_static_ch = 0;
sig_base.use_gpu = 0;
sig_base.precompute_wav = 0;
sig_base.freq_dep = 0;

% sig_direct = sigma(eps_result, sig_base, sys, options, syms);

sig_cauchy_input = sig_base;
sig_cauchy_input.isdf.enable = true;
sig_cauchy_input.isdf.algorithm = 'cauchy_cohsex';
sig_cauchy_input.isdf.sample_method = 'qrcp';
sig_cauchy_input.isdf.rank = sig_base.nbnd;
sig_cauchy_input.isdf.seed = 0;
sig_cauchy_input.isdf.cauchy_method = 'direct';
sig_cauchy_input.isdf.cauchy_froErr = 1e-8;
sig_cauchy_input.isdf.cauchy_MaxIter = 12;

sig_cauchy = sigma(eps_result, sig_cauchy_input, sys, options, syms);

sig_relative_error = norm(sig_direct.sig(:) - sig_cauchy.sig(:)) / max(1, norm(sig_direct.sig(:)));
eqp_relative_error = norm(sig_direct.eqp0(:) - sig_cauchy.eqp0(:)) / max(1, norm(sig_direct.eqp0(:)));

fprintf('AgBr spinor Cauchy COHSEX sig relative error = %.3e\n', sig_relative_error);
fprintf('AgBr spinor Cauchy COHSEX eqp0 relative error = %.3e\n', eqp_relative_error);

assert(sig_relative_error < 1e-8, ...
    'AgBr spinor Cauchy COHSEX sig validation failed: %.3e', sig_relative_error);
assert(eqp_relative_error < 1e-8, ...
    'AgBr spinor Cauchy COHSEX eqp0 validation failed: %.3e', eqp_relative_error);
