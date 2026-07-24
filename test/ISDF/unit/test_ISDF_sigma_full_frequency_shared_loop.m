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
assert(size(eps_result.inv{1}, 1) == 1, ...
    'Full-frequency shared-loop fixture must retain nmtx = 1.');

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
sig_input.exact_static_ch = 1;
sig_input.use_gpu = 0;
sig_input.precompute_wav = 0;

sig_direct = sigma(eps_result, sig_input, sys, options, syms);
local_assert_dynamic_shapes(sig_direct, sys, sig_input);

sig_matrix_input = sig_input;
sig_matrix_input.isdf.enable = true;
sig_matrix_input.isdf.algorithm = 'matrix_elements';
sig_matrix_input.isdf.sample_method = 'qrcp';
sig_matrix_input.isdf.rank = sig_input.nbnd;
sig_matrix_input.isdf.seed = 0;
sig_matrix = sigma(eps_result, sig_matrix_input, sys, options, syms);
local_assert_dynamic_shapes(sig_matrix, sys, sig_input);

fields = {'sig', 'cor', 'eqp0', 'eqp1'};
errors = zeros(size(fields));
for ifield = 1:numel(fields)
    field = fields{ifield};
    assert(isequal(size(sig_direct.(field)), size(sig_matrix.(field))), ...
        '%s output shape changed between direct and matrix-elements.', field);
    errors(ifield) = norm( ...
        sig_direct.(field)(:) - sig_matrix.(field)(:)) / ...
        max(1, norm(sig_direct.(field)(:)));
    assert(errors(ifield) < 1e-10, ...
        'Full-frequency %s relative error %.3e exceeds tolerance.', ...
        field, errors(ifield));
end
fprintf(['ISDF shared sigma full-frequency exact-CH test passed. ', ...
    'sig = %.3e, cor = %.3e, eqp0 = %.3e, eqp1 = %.3e\n'], ...
    errors(1), errors(2), errors(3), errors(4));

function local_assert_dynamic_shapes(sig, sys, input)
fields = {'sig', 'cor', 'eqp0', 'eqp1', 'neqp1'};
for ifield = 1:numel(fields)
    field = fields{ifield};
    assert(isfield(sig, field), ...
        'Full-frequency sigma result is missing %s.', field);
    assert(all(isfinite(sig.(field)(:))), ...
        'Full-frequency sigma result %s contains non-finite values.', field);
    assert(size(sig.(field), 1) == 1 && ...
        size(sig.(field), 2) == sys.nkpts && ...
        size(sig.(field), 3) == sys.nspin, ...
        'Unexpected public output shape for %s.', field);
end
assert(sig.nfreq_grid == ...
    2 * fix(input.max_freq_eval / input.delta_freq_eval) + 1);
assert(isequal(size(sig.freq_grid), [1, sig.nfreq_grid]));
end
