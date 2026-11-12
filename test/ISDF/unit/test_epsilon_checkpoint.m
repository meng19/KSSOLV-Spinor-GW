% Verify epsilon checkpoints can be written and recovered for sigma.
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(genpath(fullfile(repo_root, 'src')));

checkpoint_file = [tempname, '.mat'];
eps = struct('nmtx', 2, 'freq', 0, 'nfreq', 1);
eps.inv = {eye(2)};
request_eps = struct('cutoff', 8, 'isdf', struct('enable', false));
sys = struct('ecut', 10, 'name', 'checkpoint-test');
options = struct('kpts', [0, 0, 0]);
syms = struct('name', 'identity');

epsilon_checkpoint_save(checkpoint_file, eps, request_eps, sys, options, syms);
[loaded_eps, loaded_sys, loaded_options, loaded_syms, loaded_request] = ...
    epsilon_checkpoint_load(checkpoint_file);
assert(isequal(loaded_eps, eps));
assert(isequal(loaded_sys, sys));
assert(isequal(loaded_options, options));
assert(isequal(loaded_syms, syms));
assert(isequal(loaded_request, request_eps));

[cache_hit, cached_eps] = epsilon_checkpoint_restore( ...
    checkpoint_file, request_eps, sys, options, syms);
assert(cache_hit && isequal(cached_eps, eps));
isdf_request = request_eps;
isdf_request.isdf.enable = true;
[cache_hit, ~] = epsilon_checkpoint_restore( ...
    checkpoint_file, isdf_request, sys, options, syms);
assert(~cache_hit, 'ISDF and direct epsilon caches must not be shared.');

sigma_file = fullfile(repo_root, 'src', 'GW', 'sigma', 'sigma.m');
sigma_source = fileread(sigma_file);
assert(contains(sigma_source, 'sigma_checkpoint_inputs'));
checkpoint_inputs_file = fullfile(repo_root, 'src', 'GW', 'sigma', ...
    'sigma_checkpoint_inputs.m');
checkpoint_inputs_source = fileread(checkpoint_inputs_file);
assert(contains(checkpoint_inputs_source, ...
    'epsilon_checkpoint_load(checkpoint_file)'));
assert(contains(checkpoint_inputs_source, 'sig.epsilon_checkpoint_file'));

delete(checkpoint_file);
fprintf('Epsilon checkpoint test passed.\n');
