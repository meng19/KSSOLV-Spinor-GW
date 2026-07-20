% Validate optional ISDF sigma path for MoS2 against the existing direct path.
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

read_vxc = 1;
[sys, options, syms] = read_qe_gw('.\example\qe_data\mos2_222_spinor', read_vxc);
[sys, options] = gw_setup(sys, options);

eps_input.nbnd = 29;
eps_input.nv = options.nv;
eps_input.nc = eps_input.nbnd - eps_input.nv;
eps_input.freq_dep = 2;
eps_input.freq_dep_method = 2;
eps_input.freq_cutoff = 200;
eps_input.delta_freq = 15;
eps_input.nfreq_imag = 15;
eps_input.cutoff = 10;
eps_input.coul_cut = 'spherical_truncation';
eps_input.coul_cutoff = 10;
eps_input.use_gpu = 0;
eps_input.save_mem = 0;
eps_input.precompute_wav = 0;

eps_result = epsilon(sys, options, syms, eps_input);

sig_base.nbnd = 29;
sig_base.ndiag_min = 29;
sig_base.ndiag_max = 29;
sig_base.freq_dep = 2;
sig_base.freq_dep_method = 2;
sig_base.freq_grid_shift = 2;
sig_base.max_freq_eval = 2;
sig_base.delta_freq_eval = 0.2;
sig_base.cd_int_method = 0;
sig_base.coul_cut = 'spherical_truncation';
sig_base.coul_cutoff = 10;
sig_base.no_symmetries_q_grid = 1;
sig_base.exact_static_ch = 1;
sig_base.use_gpu = 0;
sig_base.precompute_wav = 0;

% sig_direct = sigma(eps_result, sig_base, sys, options, syms);

sig_isdf_input = sig_base;
sig_isdf_input.isdf.enable = true;
sig_isdf_input.isdf.algorithm = 'matrix_elements';
sig_isdf_input.isdf.sample_method = 'qrcp';
sig_isdf_input.isdf.rank = sig_base.nbnd;
sig_isdf_input.isdf.seed = 0;

sig_isdf = sigma(eps_result, sig_isdf_input, sys, options, syms);

sig_relative_error = norm(sig_direct.sig(:) - sig_isdf.sig(:)) / max(1, norm(sig_direct.sig(:)));
eqp_relative_error = norm(sig_direct.eqp0(:) - sig_isdf.eqp0(:)) / max(1, norm(sig_direct.eqp0(:)));

fprintf('MoS2 sigma ISDF validation sig relative error = %.3e\n', sig_relative_error);
fprintf('MoS2 sigma ISDF validation eqp0 relative error = %.3e\n', eqp_relative_error);

assert(sig_relative_error < 1e-8, ...
    'MoS2 sigma ISDF validation failed for sig: relative error %.3e', sig_relative_error);
assert(eqp_relative_error < 1e-8, ...
    'MoS2 sigma ISDF validation failed for eqp0: relative error %.3e', eqp_relative_error);
