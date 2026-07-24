clc;
clear;

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

[sys, options, syms] = read_qe_gw('.\example\qe_data\AgBr', 1);
[sys, options] = gw_setup(sys, options);

eps_input.nbnd = options.nv + 2;
eps_input.nv = options.nv;
eps_input.nc = eps_input.nbnd - eps_input.nv;
eps_input.freq_dep = 0;
eps_input.cutoff = 0.01;
eps_input.coul_cut = 'cell_box_truncation';
eps_input.use_gpu = 0;
eps_input.save_mem = 1;
eps_input.precompute_wav = 0;
eps_result = epsilon(sys, options, syms, eps_input);

% A frequency-independent screened kernel is sufficient to exercise the
% full-frequency sigma orchestration without constructing a costly dynamic
% epsilon fixture.
nfreq_real = 3;
nfreq_imag = 3;
tmp = 0:(nfreq_imag - 1);
freq_imag = -2i * 13.6056923 * tmp ./ (tmp - nfreq_imag);
eps_result.freq_dep = 2;
eps_result.freq_dep_method = 2;
eps_result.nfreq = nfreq_real + nfreq_imag;
eps_result.nfreq_imag = nfreq_imag;
eps_result.freq = [0, 15, 30, freq_imag];
for iq = 1:numel(eps_result.inv)
    eps_result.inv{iq} = repmat( ...
        eps_result.inv{iq}, 1, 1, eps_result.nfreq);
end

sig_input.nbnd = eps_input.nbnd;
sig_input.ndiag_min = options.nv;
sig_input.ndiag_max = options.nv;
sig_input.freq_dep = 2;
sig_input.freq_dep_method = 2;
sig_input.freq_grid_shift = 2;
sig_input.max_freq_eval = 2;
sig_input.delta_freq_eval = 0.2;
sig_input.cd_int_method = 0;
sig_input.coul_cut = 'cell_box_truncation';
sig_input.no_symmetries_q_grid = 0;
sig_input.exact_static_ch = 0;
sig_input.use_gpu = 0;
sig_input.precompute_wav = 0;

sig_result = sigma(eps_result, sig_input, sys, options, syms);

assert(all(isfinite(sig_result.sig(:))));
assert(all(isfinite(sig_result.eqp0(:))));
fprintf('ISDF shared sigma full-frequency dynamic test passed.\n');
