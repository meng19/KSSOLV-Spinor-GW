clc;
clear;

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(17, 'twister');

nmu = 5;
nv = 3;
nc = 4;
phi = randn(nmu, nv) + 1i * randn(nmu, nv);
psi = randn(nmu, nc) + 1i * randn(nmu, nc);
ev_occ = [-5.0; -3.0; -1.0];
ev_unocc = [0.8; 2.0; 3.5; 5.0];

direct = zeros(nmu, nmu);
for iv = 1:nv
    for ic = 1:nc
        c = conj(phi(:, iv)) .* psi(:, ic);
        direct = direct + (c * c') / (ev_occ(iv) - ev_unocc(ic));
    end
end

options.method = 'direct';
space.left_mu_components = {phi};
space.right_mu_components = {psi};
polar = isdf.polarizability(space, ev_occ, ev_unocc, options);
actual = polar.coeff;

relative_error = norm(actual - direct, 'fro') / max(1, norm(direct, 'fro'));
assert(relative_error < 1e-12, ...
    'ISDF COmegaCstar direct helper mismatch: %.3e', relative_error);

fprintf('ISDF COmegaCstar direct test passed. relative_error = %.3e\n', relative_error);

cauchy_options.method = 'cauchy';
cauchy_options.froErr = 1e-9;
cauchy_options.MaxIter = 12;
cauchy_polar = isdf.polarizability( ...
    space, ev_occ, ev_unocc, cauchy_options);
cauchy_actual = cauchy_polar.coeff;
cauchy_info = cauchy_polar.info;

cauchy_relative_error = norm(cauchy_actual - direct, 'fro') / max(1, norm(direct, 'fro'));
assert(cauchy_relative_error < 1e-6, ...
    'ISDF COmegaCstar Cauchy helper mismatch: %.3e', cauchy_relative_error);

component_polar = isdf.polarizability( ...
    space, ev_occ, ev_unocc, options);
component_actual = component_polar.coeff;
wrapper_relative_error = norm(actual - component_actual, 'fro') / max(1, norm(component_actual, 'fro'));
assert(wrapper_relative_error < 1e-12, ...
    'Scalar COmegaCstar should match single-component cell path: %.3e', wrapper_relative_error);

fprintf('ISDF COmegaCstar Cauchy test passed. relative_error = %.3e, iterations = %d\n', ...
    cauchy_relative_error, cauchy_info.iterations);

dynamic_freq = [0, 0.2, 0.2i, 0.5i];
dynamic_direct = isdf.polarizability(space, ev_occ, ev_unocc, ...
    struct('method', 'direct', 'freq', dynamic_freq));
dynamic_cauchy = isdf.polarizability(space, ev_occ, ev_unocc, ...
    struct('method', 'cauchy', 'freq', dynamic_freq, ...
    'froErr', 1e-9, 'MaxIter', 12));
dynamic_relative_error = norm( ...
    dynamic_cauchy.coeff(:) - dynamic_direct.coeff(:)) / ...
    max(1, norm(dynamic_direct.coeff(:)));
assert(dynamic_relative_error < 1e-6, ...
    'Dynamic Cauchy COmegaCstar mismatch: %.3e', ...
    dynamic_relative_error);
fprintf('ISDF COmegaCstar dynamic Cauchy test passed. relative_error = %.3e\n', ...
    dynamic_relative_error);

fallback_freq = 2.7;
fallback_direct = isdf.polarizability(space, ev_occ, ev_unocc, ...
    struct('method', 'direct', 'freq', fallback_freq));
fallback_cauchy = isdf.polarizability(space, ev_occ, ev_unocc, ...
    struct('method', 'cauchy', 'freq', fallback_freq, ...
    'froErr', 1e-9, 'MaxIter', 12));
fallback_relative_error = norm( ...
    fallback_cauchy.coeff(:) - fallback_direct.coeff(:)) / ...
    max(1, norm(fallback_direct.coeff(:)));
assert(fallback_relative_error < 1e-12, ...
    'Cauchy direct fallback mismatch: %.3e', ...
    fallback_relative_error);
assert(fallback_cauchy.info.fallback_count == 1);
assert(isequal(fallback_cauchy.info.fallback_pages, true));
fprintf('ISDF COmegaCstar fallback info test passed. fallback_count = %d\n', ...
    fallback_cauchy.info.fallback_count);
