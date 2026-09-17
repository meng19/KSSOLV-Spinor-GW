script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(29, 'twister');
ngrid = 14;
nmu_vc = 3;
nmu_t = 4;

zeta_g = randn(ngrid, nmu_vc) + 1i * randn(ngrid, nmu_vc);
epsilon_vcoul = 0.5 + rand(ngrid, 1);
k_mu = randn(nmu_vc, nmu_vc) + 1i * randn(nmu_vc, nmu_vc);
screened = struct('zeta_g', zeta_g, 'epsilon_vcoul', epsilon_vcoul, ...
    'k_mu', k_mu);
target = randn(ngrid, nmu_t) + 1i * randn(ngrid, nmu_t);
contract_vcoul = 0.25 + rand(ngrid, 1);

isdf.screened_kernel_cache('reset');
isdf.screened_kernel_cache('limit', 1e9);
info = isdf.screened_kernel_cache('stats');
assert(info.limit == 1e9 && info.stored == 0 && info.hits == 0 && ...
    info.misses == 0 && info.bytes == 0);

% An uncached call still matches the explicit projection and must leave the
% store untouched.
reference = isdf.screened_kernel(screened, target, contract_vcoul);
left_projector = target.' * (epsilon_vcoul .* zeta_g);
reference_direct = left_projector * k_mu * ...
    (zeta_g' * (contract_vcoul .* conj(target)));
assert(norm(reference - reference_direct, 'fro') / ...
    norm(reference_direct, 'fro') < 1e-12);
info = isdf.screened_kernel_cache('stats');
assert(info.stored == 0, 'Uncached calls must not populate the store.');

[missing, hit] = isdf.screened_kernel_cache('get', 'not-stored');
assert(isempty(missing) && ~hit);

key = 'target|nn-space-k1-q1-s1-b29|q1|n14|v3|m3|p1|g0';
kernel = isdf.screened_kernel(screened, target, contract_vcoul, key);
assert(isequal(kernel, reference), ...
    'A cached call must return the same kernel as the direct projection.');
info = isdf.screened_kernel_cache('stats');
assert(info.stored == 1 && info.misses == 2 && info.hits == 0);
assert(info.bytes == numel(kernel) * 16, ...
    'Stored complex kernels must count 16 bytes per element.');

reused = isdf.screened_kernel(screened, target, contract_vcoul, key);
assert(isequal(reused, kernel));
info = isdf.screened_kernel_cache('stats');
assert(info.hits == 1 && info.stored == 1);

% A different key recomputes the projection instead of reusing the entry.
other_key = 'target|nn-space-k1-q1-s1-b30|q1|n14|v3|m3|p1|g0';
other = isdf.screened_kernel(screened, target, contract_vcoul, other_key);
assert(isequal(other, reference));
info = isdf.screened_kernel_cache('stats');
assert(info.stored == 2 && info.misses == 3);

% The store trusts the key: sigma resets it at the start of every run, so a
% stale key can only appear when a caller reuses keys across changed data.
changed = screened;
changed.k_mu = 2 * k_mu;
stale = isdf.screened_kernel(changed, target, contract_vcoul, key);
assert(isequal(stale, kernel));

% A budget that cannot hold the kernel skips storing it instead of evicting
% or exceeding the limit.
isdf.screened_kernel_cache('reset');
isdf.screened_kernel_cache('limit', numel(kernel) * 16 - 1);
third = isdf.screened_kernel(screened, target, contract_vcoul, key);
assert(isequal(third, reference));
info = isdf.screened_kernel_cache('stats');
assert(info.stored == 0 && info.skipped == 1 && info.bytes == 0);
[~, hit] = isdf.screened_kernel_cache('get', key);
assert(~hit);

% Frequency-page (3-D k_mu) kernels follow the same protocol.
screened_pages = screened;
screened_pages.k_mu = cat(3, k_mu, 0.5 * k_mu);
page_reference = isdf.screened_kernel(screened_pages, target, contract_vcoul);
isdf.screened_kernel_cache('reset');
isdf.screened_kernel_cache('limit', 1e9);
page_kernel = isdf.screened_kernel(screened_pages, target, ...
    contract_vcoul, 'target|pages');
assert(isequal(page_kernel, page_reference));
assert(isequal(size(page_kernel), [nmu_t, nmu_t, 2]));
page_reused = isdf.screened_kernel(screened_pages, target, ...
    contract_vcoul, 'target|pages');
assert(isequal(page_reused, page_kernel));
info = isdf.screened_kernel_cache('stats');
assert(info.stored == 1 && info.hits == 1 && info.misses == 1);

% The target-free (full matrix) projection is cached with its own key.
full_reference = isdf.screened_kernel(screened, [], contract_vcoul);
full_kernel = isdf.screened_kernel(screened, [], contract_vcoul, 'full|q1');
assert(isequal(full_kernel, full_reference));
full_reused = isdf.screened_kernel(screened, [], contract_vcoul, 'full|q1');
assert(isequal(full_reused, full_kernel));
info = isdf.screened_kernel_cache('stats');
assert(info.stored == 2 && info.hits == 2 && info.misses == 2);

% Replacing an existing entry keeps a single stored kernel.
isdf.screened_kernel_cache('reset');
isdf.screened_kernel_cache('limit', 1e9);
isdf.screened_kernel_cache('put', 'duplicate', zeros(2));
isdf.screened_kernel_cache('put', 'duplicate', zeros(3));
info = isdf.screened_kernel_cache('stats');
assert(info.stored == 1 && info.bytes == 9 * 8);

% A zero budget disables storage, and invalid input is rejected.
isdf.screened_kernel_cache('limit', 0);
info = isdf.screened_kernel_cache('stats');
assert(info.limit == 0);
isdf.screened_kernel_cache('reset');

limit_error = '';
try
    isdf.screened_kernel_cache('limit', -1);
catch ME
    limit_error = ME.identifier;
end
assert(strcmp(limit_error, 'ISDF:ScreenedKernelCacheLimit'), ...
    'Expected an ISDF:ScreenedKernelCacheLimit error, got "%s".', limit_error);

argument_error = '';
try
    isdf.screened_kernel_cache('limit');
catch ME
    argument_error = ME.identifier;
end
assert(strcmp(argument_error, 'ISDF:ScreenedKernelCacheArgument'), ...
    'Expected an ISDF:ScreenedKernelCacheArgument error, got "%s".', ...
    argument_error);

action_error = '';
try
    isdf.screened_kernel_cache('unknown-action');
catch ME
    action_error = ME.identifier;
end
assert(strcmp(action_error, 'ISDF:UnknownScreenedKernelCacheAction'), ...
    'Expected an ISDF:UnknownScreenedKernelCacheAction error, got "%s".', ...
    action_error);

% Sigma must clear and configure the store for every run.
sigma_source = fileread(fullfile(repo_root, 'src', 'GW', 'sigma', 'sigma.m'));
assert(contains(sigma_source, 'isdf.screened_kernel_cache(''reset'');'));
assert(contains(sigma_source, 'isdf.screened_kernel_cache(''limit'', ...'));

fprintf('ISDF screened-kernel cache test passed.\n');
