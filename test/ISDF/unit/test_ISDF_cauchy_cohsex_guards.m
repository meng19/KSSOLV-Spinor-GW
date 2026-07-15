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
assert(strcmp(sig.isdf.algorithm, 'matrix_elements'));
assert(isfield(sig.isdf, 'cauchy_method'));
assert(strcmp(sig.isdf.cauchy_method, 'cauchy'));
assert(isfield(sig.isdf, 'cauchy_froErr'));
assert(isfield(sig.isdf, 'cauchy_MaxIter'));

fprintf('ISDF Cauchy COHSEX default guard test passed.\n');
