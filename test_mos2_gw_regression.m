% Compare accepted MoS2 multi-k GW calculations against Markdown baselines.
%
% Optional variables before run():
%   mos2_cases = {'mos2_222_cohsex'};  % or {'mos2_444_fullfreq_gw'}
%   mos2_use_gpu = false;              % override eps/sig use_gpu
%   mos2_abs_tol = 1e-7;
%   mos2_rel_tol = 1e-10;

if ~exist('mos2_cases', 'var')
    mos2_cases = {'mos2_222_cohsex', 'mos2_444_fullfreq_gw'};
end
if ~exist('mos2_use_gpu', 'var')
    mos2_use_gpu = [];
end
if ~exist('mos2_abs_tol', 'var')
    mos2_abs_tol = 1e-7;
end
if ~exist('mos2_rel_tol', 'var')
    mos2_rel_tol = 1e-10;
end

script_dir = fileparts(mfilename('fullpath'));
repo_root = script_dir;
baseline_file = fullfile(repo_root, 'docs', 'mos2-gw-baselines.md');

old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
addpath(repo_root);

randn('state', 0);
rand('state', 0);
KSSOLV_startup;

for icase = 1:numel(mos2_cases)
    case_id = mos2_cases{icase};
    baseline = load_mos2_baseline(baseline_file, case_id);
    actual = run_mos2_case(case_id, mos2_use_gpu);
    compare_mos2_result(case_id, actual, baseline, mos2_abs_tol, ...
        mos2_rel_tol);
end

fprintf('MoS2 GW regression baselines passed for %d case(s).\n', ...
    numel(mos2_cases));

function result = run_mos2_case(case_id, use_gpu_override)
switch case_id
    case 'mos2_222_cohsex'
        result = run_mos2_222_cohsex(use_gpu_override);
    case 'mos2_444_fullfreq_gw'
        result = run_mos2_444_fullfreq_gw(use_gpu_override);
    otherwise
        error('Unknown MoS2 regression case: %s', case_id);
end
end

function result = run_mos2_222_cohsex(use_gpu_override)
[sys, options, syms] = read_qe_gw('.\example\qe_data\mos2_222_spinor', 1);
[sys, options] = gw_setup(sys, options);

eps_input.nbnd = 29;
eps_input.nv = options.nv;
eps_input.nc = eps_input.nbnd - eps_input.nv;
eps_input.freq_dep = 0;
eps_input.cutoff = 10;
eps_input.coul_cut = 'spherical_truncation';
eps_input.coul_cutoff = 10;
eps_input.use_gpu = 0;
eps_input.save_mem = 0;
eps_input.precompute_wav = 0;
eps_input = apply_gpu_override(eps_input, use_gpu_override);

eps_result = epsilon(sys, options, syms, eps_input);

sig_input.nbnd = 29;
sig_input.ndiag_min = 26;
sig_input.ndiag_max = 29;
sig_input.freq_dep = 0;
sig_input.coul_cut = 'spherical_truncation';
sig_input.coul_cutoff = 10;
sig_input.no_symmetries_q_grid = 1;
sig_input.exact_static_ch = 1;
sig_input.use_gpu = 0;
sig_input.precompute_wav = 0;
sig_input = apply_gpu_override(sig_input, use_gpu_override);

sig_result = sigma(eps_result, sig_input, sys, options, syms);
result = collect_mos2_result('mos2_222_cohsex', sig_result, sys);
end

function result = run_mos2_444_fullfreq_gw(use_gpu_override)
[sys, options, syms] = read_qe_gw( ...
    '.\example\qe_data\mos2_444_spinor_fullfreq', 1);
[sys, options] = gw_setup(sys, options);

eps_input.nbnd = 29;
eps_input.nv = options.nv;
eps_input.nc = eps_input.nbnd - eps_input.nv;
eps_input.freq_dep = 2;
eps_input.freq_dep_method = 2;
eps_input.freq_cutoff = 200;
eps_input.delta_freq = 100;
eps_input.nfreq_imag = 2;
eps_input.cutoff = 5;
eps_input.coul_cut = 'spherical_truncation';
eps_input.coul_cutoff = 5;
eps_input.use_gpu = 1;
eps_input.save_mem = 1;
eps_input.precompute_wav = 0;
eps_input = apply_gpu_override(eps_input, use_gpu_override);

eps_result = epsilon(sys, options, syms, eps_input);

sig_input.nbnd = 29;
sig_input.ndiag_min = 28;
sig_input.ndiag_max = 29;
sig_input.freq_dep = 2;
sig_input.freq_dep_method = 2;
sig_input.freq_grid_shift = 2;
sig_input.max_freq_eval = 2;
sig_input.delta_freq_eval = 0.2;
sig_input.cd_int_method = 0;
sig_input.coul_cut = 'spherical_truncation';
sig_input.coul_cutoff = 5;
sig_input.no_symmetries_q_grid = 1;
sig_input.exact_static_ch = 1;
sig_input.use_gpu = 1;
sig_input.precompute_wav = 0;
sig_input = apply_gpu_override(sig_input, use_gpu_override);

sig_result = sigma(eps_result, sig_input, sys, options, syms);
result = collect_mos2_result('mos2_444_fullfreq_gw', sig_result, sys);
end

function input = apply_gpu_override(input, use_gpu_override)
if ~isempty(use_gpu_override)
    input.use_gpu = double(logical(use_gpu_override));
end
end

function result = collect_mos2_result(case_id, sig_result, sys)
result.case_id = case_id;
result.sig = sig_result.sig;
result.eqp0 = sig_result.eqp0;
result.vxc = sig_result.vxc;
result.kpts = sys.kpts;
end

function baseline = load_mos2_baseline(filename, case_id)
text = fileread(filename);
marker = ['<!-- baseline:' case_id ' -->'];
start_idx = strfind(text, marker);
assert(~isempty(start_idx), 'Missing baseline marker for %s.', case_id);

remaining = text(start_idx(1):end);
tokens = regexp(remaining, '```json\s*([\s\S]*?)\s*```', 'tokens', ...
    'once');
assert(~isempty(tokens), 'Missing JSON baseline block for %s.', case_id);

baseline = jsondecode(tokens{1});
assert(strcmp(baseline.case_id, case_id), ...
    'Baseline case id mismatch: expected %s, found %s.', ...
    case_id, baseline.case_id);
end

function compare_mos2_result(case_id, actual, expected, abs_tol, rel_tol)
fields = {'sig', 'eqp0', 'vxc', 'kpts'};
for ifield = 1:numel(fields)
    field = fields{ifield};
    actual_value = actual.(field);
    expected_value = expected.(field);
    assert(isequal(size(actual_value), size(expected_value)), ...
        '%s.%s size mismatch: actual %s, expected %s.', case_id, field, ...
        mat2str(size(actual_value)), mat2str(size(expected_value)));

    diff_value = actual_value - expected_value;
    max_abs_error = max(abs(diff_value(:)));
    rel_error = norm(diff_value(:)) / max(1, norm(expected_value(:)));
    fprintf('%s %-5s max_abs = %.3e, rel = %.3e\n', case_id, field, ...
        max_abs_error, rel_error);
    assert(max_abs_error <= abs_tol || rel_error <= rel_tol, ...
        '%s.%s differs from baseline: max_abs %.3e, rel %.3e.', ...
        case_id, field, max_abs_error, rel_error);
end
end
