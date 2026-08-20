% Validate full-frequency MoS2 2x2x2 spinor reduced-basis ISDF paths.
clc;
clear;
close all;
rng(0, 'twister');

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

[sys, options, syms] = read_qe_gw( ...
    '.\example\qe_data\mos2_222_spinor', 1);
[sys, options] = gw_setup(sys, options);
assert(sys.nkpts > 1, 'This validation requires a multi-k data set.');
assert(sys.nspinor > 1, 'This validation requires spinor wavefunctions.');

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

tic
eps_direct = epsilon(sys, options, syms, eps_input);

sig_input.nbnd = eps_input.nbnd;
sig_input.ndiag_min = eps_input.nbnd;
sig_input.ndiag_max = eps_input.nbnd;
sig_input.freq_dep = 2;
sig_input.freq_dep_method = 2;
sig_input.freq_grid_shift = 2;
sig_input.max_freq_eval = 2;
sig_input.delta_freq_eval = 0.2;
sig_input.cd_int_method = 0;
sig_input.coul_cut = 'spherical_truncation';
sig_input.coul_cutoff = 10;
sig_input.no_symmetries_q_grid = 1;
sig_input.exact_static_ch = 1;
sig_input.use_gpu = 0;
sig_input.precompute_wav = 0;

sig_direct = sigma(eps_direct, sig_input, sys, options, syms);
toc

eps_reduced_input = eps_input;
eps_reduced_input.isdf.enable = true;
eps_reduced_input.isdf.algorithm = 'reduced_basis';
eps_reduced_input.isdf.output = 'screened_w';
% eps_reduced_input.isdf.output = 'both';
eps_reduced_input.isdf.sample_method = 'qrcp';
eps_reduced_input.isdf.rank = eps_input.nv * eps_input.nc;
eps_reduced_input.isdf.reduced_solver = 'direct';
eps_reduced_input.isdf.seed = 0;

tic
eps_reduced = epsilon(sys, options, syms, eps_reduced_input);

% eps_error = 0;
% for iq = 1:sys.nkpts
%     eps_error = max(eps_error, local_relative_error( ...
%         eps_direct.inv{iq}, eps_reduced.inv{iq}));
%     assert(~isempty(eps_reduced.isdf_screened_w{iq}), ...
%         'Reduced screened W is missing at q-point %d.', iq);
% end
% fprintf(['MoS2 2x2x2 spinor full-frequency reduced epsilon ', ...
%     'relative error = %.3e\n'], eps_error);
% assert(eps_error < 1e-8, ...
%     'Full-frequency reduced epsilon validation failed: %.3e', ...
%     eps_error);

sig_reduced_input = sig_input;
sig_reduced_input.isdf.enable = true;
sig_reduced_input.isdf.algorithm = 'reduced_basis';
sig_reduced_input.isdf.sample_method = 'qrcp';
sig_reduced_input.isdf.rank = sig_input.nbnd;
sig_reduced_input.isdf.seed = 0;

sig_reduced = sigma(eps_reduced, sig_reduced_input, sys, options, syms);
toc

sig_error = local_relative_error(sig_direct.sig, sig_reduced.sig);
eqp0_error = local_relative_error(sig_direct.eqp0, sig_reduced.eqp0);
eqp1_error = local_relative_error(sig_direct.eqp1, sig_reduced.eqp1);
fprintf(['MoS2 2x2x2 spinor full-frequency reduced sigma ', ...
    'relative error = %.3e\n'], sig_error);
fprintf(['MoS2 2x2x2 spinor full-frequency reduced eqp0 ', ...
    'relative error = %.3e\n'], eqp0_error);
fprintf(['MoS2 2x2x2 spinor full-frequency reduced eqp1 ', ...
    'relative error = %.3e\n'], eqp1_error);
assert(sig_error < 1e-8, ...
    'Full-frequency reduced sigma validation failed: %.3e', ...
    sig_error);
assert(eqp0_error < 1e-8, ...
    'Full-frequency reduced eqp0 validation failed: %.3e', ...
    eqp0_error);
assert(eqp1_error < 1e-8, ...
    'Full-frequency reduced eqp1 validation failed: %.3e', ...
    eqp1_error);

% ---- Relative error helper ----

function error_value = local_relative_error(reference, candidate)
numerator = norm(reference(:) - candidate(:));
denominator = max(1, norm(reference(:)));
error_value = numerator / denominator;
end
