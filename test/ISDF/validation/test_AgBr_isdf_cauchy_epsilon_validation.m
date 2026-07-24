% Validate component ISDF reduced-basis epsilon output.
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
[sys, options, syms] = read_qe_gw('.\example\qe_data\AgBr', read_vxc);
[sys, options] = gw_setup(sys, options);

eps_base.nbnd = options.nv + 2;
eps_base.nv = options.nv;
eps_base.nc = eps_base.nbnd - eps_base.nv;
eps_base.cutoff = 2;
eps_base.coul_cut = 'cell_box_truncation';
eps_base.use_gpu = 0;
eps_base.save_mem = 1;
eps_base.precompute_wav = 0;
eps_base.freq_dep = 0;

eps_direct = epsilon(sys, options, syms, eps_base);

eps_cauchy_input = eps_base;
eps_cauchy_input.isdf.enable = true;
eps_cauchy_input.isdf.algorithm = 'reduced_basis';
eps_cauchy_input.isdf.sample_method = 'qrcp';
eps_cauchy_input.isdf.rank = 8;
eps_cauchy_input.isdf.seed = 0;
eps_cauchy_input.isdf.reduced_solver = 'cauchy';
eps_cauchy_input.isdf.cauchy_froErr = 1e-8;
eps_cauchy_input.isdf.cauchy_MaxIter = 12;

eps_cauchy = epsilon(sys, options, syms, eps_cauchy_input);

assert(~isfield(eps_cauchy, 'inv'), ...
    'Reduced-basis epsilon should not store a full G-space inverse.');
assert(isfield(eps_cauchy, 'isdf_screened_w'), ...
    'Reduced-basis epsilon must store reduced screened W.');
assert(isfield(eps_cauchy.isdf_screened_w{1, 1}, 'k_mu'), ...
    'Reduced screened W must contain the reduced-space K matrix.');
assert(strcmp(eps_cauchy.isdf_reduced_info{1, 1}.method, 'cauchy'));

eps_both_input = eps_cauchy_input;
eps_both_input.isdf.output = 'both';
eps_both = epsilon(sys, options, syms, eps_both_input);
assert(isfield(eps_both, 'inv'), ...
    'Reduced-basis epsilon output="both" must store full eps.inv.');
assert(isfield(eps_both, 'isdf_screened_w'), ...
    'Reduced-basis epsilon output="both" must store reduced screened W.');

max_relative_error = 0;
for iq = 1:length(eps_direct.inv)
    numerator = norm(eps_direct.inv{iq}(:) - eps_both.inv{iq}(:));
    denominator = max(1, norm(eps_direct.inv{iq}(:)));
    max_relative_error = max(max_relative_error, numerator / denominator);
end

fprintf('AgBr component reduced-basis epsilon full-inverse relative error = %.3e\n', max_relative_error);
assert(max_relative_error < 1e-2, ...
    'AgBr reduced-basis full-inverse validation failed: relative error %.3e', max_relative_error);

fprintf('AgBr component reduced-basis epsilon validation passed.\n');
