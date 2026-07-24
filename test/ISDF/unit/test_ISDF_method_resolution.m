script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

assert(strcmp(gw_resolve_method([], 'epsilon'), 'direct'));
assert(strcmp(gw_resolve_method(struct(), 'epsilon'), 'direct'));
assert(strcmp(gw_resolve_method('legacy-value', 'sigma'), 'direct'));
assert(strcmp(gw_resolve_method(struct( ...
    'enable', false, 'algorithm', 'invalid', 'output', 17), ...
    'epsilon'), 'direct'));
assert(strcmp(gw_resolve_method(struct( ...
    'enable', false, 'algorithm', 'invalid', 'reduced_solver', []), ...
    'sigma'), 'direct'));
assert(strcmp(gw_resolve_method(struct( ...
    'enable', true, 'algorithm', 'matrix_elements'), 'epsilon'), ...
    'matrix_elements'));
assert(strcmp(gw_resolve_method(struct( ...
    'enable', true, 'algorithm', 'REDUCED_BASIS'), 'sigma'), ...
    'reduced_basis'));

try
    gw_resolve_method(struct('enable', true, 'algorithm', 'bad'), 'epsilon');
    error('TEST:ExpectedFailure', 'Unknown epsilon method did not fail.');
catch ME
    assert(strcmp(ME.identifier, 'ISDF:UnknownEpsilonAlgorithm'));
end

try
    gw_resolve_method(struct('enable', true, 'algorithm', 'bad'), 'sigma');
    error('TEST:ExpectedFailure', 'Unknown sigma method did not fail.');
catch ME
    assert(strcmp(ME.identifier, 'ISDF:UnknownSigmaAlgorithm'));
end

fprintf('ISDF method resolution test passed.\n');
