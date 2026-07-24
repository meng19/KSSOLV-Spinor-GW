script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
epsilon_file = fullfile(repo_root, 'src', 'GW', 'epsilon', 'epsilon.m');
source = fileread(epsilon_file);

assert(contains(source, 'ops = epsilon_ops(ctx);'));
assert(contains(source, 'contribution = ops.evaluate(block);'));
assert(contains(source, 'acc = ops.accumulate(acc, contribution, block);'));
assert(~contains(source, 'use_isdf = isfield'), ...
    'epsilon.m must not resolve ISDF mode inside its main loop.');
assert(~contains(source, 'isdf_epsilon_reduced_basis'));
assert(isempty(regexp(source, ...
    'isdf_epsilon_reduced_basis\([^;]+;\s*return', 'once')));
assert(exist(fullfile(repo_root, 'src', 'GW', 'ISDF', ...
    'isdf_epsilon_reduced_basis.m'), 'file') == 0);

accumulate_position = strfind(source, ...
    'acc = ops.accumulate(acc, contribution, block);');
progress_position = strfind(source, ...
    'print_progress(current_bands_for_k, total_bands_for_k, ...');
assert(accumulate_position < progress_position, ...
    'Epsilon progress must report completed block work.');

fprintf('ISDF shared epsilon loop structure test passed.\n');
