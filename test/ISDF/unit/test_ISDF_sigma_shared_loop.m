script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
sigma_file = fullfile(repo_root, 'src', 'GW', 'sigma', 'sigma.m');
source = fileread(sigma_file);

assert(contains(source, 'ops = sigma_ops(ctx);'));
assert(contains(source, 'matrix_elements = ops.matrix_elements(block);'));
assert(contains(source, 'contribution = ops.contract(block, matrix_elements);'));
assert(isempty(regexp(source, ...
    'if\s+use_isdf_matrix_elements', 'once')), ...
    'sigma.m must not select matrix construction inside the band loop.');
assert(~contains(source, 'isdf_sigma_reduced_basis'));
assert(isempty(regexp(source, ...
    'isdf_sigma_reduced_basis\([^;]+;\s*return', 'once')));
assert(exist(fullfile(repo_root, 'src', 'GW', 'ISDF', ...
    'isdf_sigma_reduced_basis.m'), 'file') == 0);

fprintf('ISDF shared sigma loop structure test passed.\n');
