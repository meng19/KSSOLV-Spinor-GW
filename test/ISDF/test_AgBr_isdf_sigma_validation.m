% Validate optional ISDF sigma path against the existing direct path.
clc;
clear all;
close all;
randn('state', 0);
rand('state', 0);

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
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

eps_result = epsilon(sys, options, syms, eps_input);

sig_base.nbnd = eps_input.nbnd;
sig_base.ndiag_min = options.nv;
sig_base.ndiag_max = options.nv;
sig_base.coul_cut = 'cell_box_truncation';
sig_base.no_symmetries_q_grid = 0;
sig_base.exact_static_ch = 0;
sig_base.use_gpu = 0;
sig_base.precompute_wav = 0;

sig_direct = sigma(eps_result, sig_base, sys, options, syms);

sig_isdf_input = sig_base;
sig_isdf_input.isdf.enable = true;
sig_isdf_input.isdf.sample_method = 'qrcp';
sig_isdf_input.isdf.rank = sig_base.nbnd;
sig_isdf_input.isdf.seed = 0;

sig_isdf = sigma(eps_result, sig_isdf_input, sys, options, syms);

sig_relative_error = norm(sig_direct.sig(:) - sig_isdf.sig(:)) / max(1, norm(sig_direct.sig(:)));
eqp_relative_error = norm(sig_direct.eqp0(:) - sig_isdf.eqp0(:)) / max(1, norm(sig_direct.eqp0(:)));

fprintf('AgBr sigma ISDF validation sig relative error = %.3e\n', sig_relative_error);
fprintf('AgBr sigma ISDF validation eqp0 relative error = %.3e\n', eqp_relative_error);

assert(sig_relative_error < 1e-8, ...
    'AgBr sigma ISDF validation failed for sig: relative error %.3e', sig_relative_error);
assert(eqp_relative_error < 1e-8, ...
    'AgBr sigma ISDF validation failed for eqp0: relative error %.3e', eqp_relative_error);
