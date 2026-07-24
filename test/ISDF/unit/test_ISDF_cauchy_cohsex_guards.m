script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

sig = struct();
sig = sigma_set_defaults(sig);
assert(isfield(sig, 'isdf'));
assert(isfield(sig.isdf, 'algorithm'));
assert(strcmp(sig.isdf.algorithm, 'reduced_basis'));
assert(isfield(sig.isdf, 'reduced_solver'));
assert(strcmp(sig.isdf.reduced_solver, 'cauchy'));
assert(isfield(sig.isdf, 'cauchy_froErr'));
assert(isfield(sig.isdf, 'cauchy_MaxIter'));

sig_matrix = struct();
sig_matrix.isdf.algorithm = 'matrix_elements';
sig_matrix = sigma_set_defaults(sig_matrix);
assert(strcmp(sig_matrix.isdf.algorithm, 'matrix_elements'));

sig_fullfreq = struct();
sig_fullfreq.freq_dep = 2;
sig_fullfreq = sigma_set_defaults(sig_fullfreq);
assert(strcmp(sig_fullfreq.isdf.algorithm, 'matrix_elements'));

eps = struct();
eps = epsilon_set_defaults(eps);
assert(isfield(eps, 'isdf'));
assert(isfield(eps.isdf, 'algorithm'));
assert(strcmp(eps.isdf.algorithm, 'reduced_basis'));
assert(isfield(eps.isdf, 'reduced_solver'));
assert(strcmp(eps.isdf.reduced_solver, 'cauchy'));
assert(isfield(eps.isdf, 'output'));
assert(strcmp(eps.isdf.output, 'screened_w'));
assert(isfield(eps.isdf, 'cauchy_froErr'));
assert(isfield(eps.isdf, 'cauchy_MaxIter'));
assert(~isfield(eps.isdf, 'store_full_inverse'));

eps_matrix = struct();
eps_matrix.isdf.algorithm = 'matrix_elements';
eps_matrix = epsilon_set_defaults(eps_matrix);
assert(strcmp(eps_matrix.isdf.algorithm, 'matrix_elements'));

eps_fullfreq = struct();
eps_fullfreq.freq_dep = 2;
eps_fullfreq = epsilon_set_defaults(eps_fullfreq);
assert(strcmp(eps_fullfreq.isdf.algorithm, 'matrix_elements'));

fprintf('ISDF reduced default guard test passed.\n');
