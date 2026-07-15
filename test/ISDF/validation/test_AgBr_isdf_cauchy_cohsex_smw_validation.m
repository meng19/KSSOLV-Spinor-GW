% Validate spinor ISDF Cauchy COHSEX with SMW-only epsilon storage.
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

eps_direct = epsilon(sys, options, syms, eps_input);

sig_base.nbnd = eps_input.nbnd;
sig_base.ndiag_min = options.nv;
sig_base.ndiag_max = options.nv;
sig_base.coul_cut = 'cell_box_truncation';
sig_base.no_symmetries_q_grid = 0;
sig_base.exact_static_ch = 0;
sig_base.use_gpu = 0;
sig_base.precompute_wav = 0;
sig_base.freq_dep = 0;

sig_direct = sigma(eps_direct, sig_base, sys, options, syms);

sig_smw_input = sig_base;
sig_smw_input.isdf.enable = true;
sig_smw_input.isdf.algorithm = 'cauchy_cohsex';
sig_smw_input.isdf.sample_method = 'qrcp';
sig_smw_input.isdf.rank = sig_base.nbnd;
sig_smw_input.isdf.seed = 0;
sig_smw_input.isdf.cauchy_method = 'cauchy';
sig_smw_input.isdf.cauchy_froErr = 1e-8;
sig_smw_input.isdf.cauchy_MaxIter = 12;

eps_cauchy_input = eps_input;
eps_cauchy_input.isdf.enable = true;
eps_cauchy_input.isdf.algorithm = 'cauchy_polarizability';
eps_cauchy_input.isdf.sample_method = 'qrcp';
eps_cauchy_input.isdf.rank = eps_input.nv * eps_input.nc;
eps_cauchy_input.isdf.seed = 0;
eps_cauchy_input.isdf.cauchy_method = 'cauchy';
eps_cauchy_input.isdf.cauchy_froErr = 1e-8;
eps_cauchy_input.isdf.cauchy_MaxIter = 12;
eps_cauchy_input.isdf.store_full_inverse = true;

eps_cauchy_full = epsilon(sys, options, syms, eps_cauchy_input);
eps_cauchy_full_inverse_only = rmfield(eps_cauchy_full, 'isdf_screened');

sig_cauchy_full = sigma(eps_cauchy_full_inverse_only, sig_smw_input, sys, options, syms);

eps_smw_input = eps_cauchy_input;
eps_smw_input.isdf.store_full_inverse = false;

eps_smw = epsilon(sys, options, syms, eps_smw_input);
assert(~isfield(eps_smw, 'inv'), 'SMW-only epsilon should not store full eps.inv.');
assert(isfield(eps_smw, 'isdf_screened'), 'SMW-only epsilon must store reduced screened object.');

sig_smw = sigma(eps_smw, sig_smw_input, sys, options, syms);

direct_sig_relative_error = norm(sig_direct.sig(:) - sig_smw.sig(:)) / max(1, norm(sig_direct.sig(:)));
direct_eqp_relative_error = norm(sig_direct.eqp0(:) - sig_smw.eqp0(:)) / max(1, norm(sig_direct.eqp0(:)));
sig_relative_error = norm(sig_cauchy_full.sig(:) - sig_smw.sig(:)) / max(1, norm(sig_cauchy_full.sig(:)));
eqp_relative_error = norm(sig_cauchy_full.eqp0(:) - sig_smw.eqp0(:)) / max(1, norm(sig_cauchy_full.eqp0(:)));

fprintf('AgBr spinor Cauchy COHSEX SMW sig relative error = %.3e\n', sig_relative_error);
fprintf('AgBr spinor Cauchy COHSEX SMW eqp0 relative error = %.3e\n', eqp_relative_error);
fprintf('AgBr spinor Cauchy COHSEX SMW direct-reference sig relative error = %.3e\n', direct_sig_relative_error);
fprintf('AgBr spinor Cauchy COHSEX SMW direct-reference eqp0 relative error = %.3e\n', direct_eqp_relative_error);

assert(sig_relative_error < 1e-8, ...
    'AgBr spinor Cauchy COHSEX SMW sig validation failed: %.3e', sig_relative_error);
assert(eqp_relative_error < 1e-8, ...
    'AgBr spinor Cauchy COHSEX SMW eqp0 validation failed: %.3e', eqp_relative_error);
