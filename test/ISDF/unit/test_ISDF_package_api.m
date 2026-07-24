script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(61, 'twister');
fftgrid = [4, 3, 2];
ngrid = prod(fftgrid);
idx_q = [1; 2; 5; 9; 13; 24];
left = {randn(ngrid, 2) + 1i * randn(ngrid, 2), ...
    randn(ngrid, 2) + 1i * randn(ngrid, 2)};
right = {randn(ngrid, 3) + 1i * randn(ngrid, 3), ...
    randn(ngrid, 3) + 1i * randn(ngrid, 3)};
options = struct('rank', 6, 'sample_method', 'qrcp', 'seed', 0);

space = isdf.build_space(left, right, idx_q, fftgrid, options);
gme = isdf.matrix_elements(left, right, idx_q, fftgrid, options);
reconstructed = space.zeta_g * space.product_mu;
assert(max(abs(gme(:) - reconstructed(:))) < 1e-12);

legacy_space = isdf_build_space(left, right, idx_q, fftgrid, options);
assert(max(abs(space.zeta_g(:) - legacy_space.zeta_g(:))) < 1e-12);
assert(isequal(space.ind_mu, legacy_space.ind_mu));

ev_occ = [-0.8; -0.2];
ev_unocc = [0.3; 0.9; 1.4];
solver = struct('method', 'direct');
polar = isdf.polarizability(space, ev_occ, ev_unocc, solver);
legacy_polar = isdf_reduced_polarizability( ...
    legacy_space, ev_occ, ev_unocc, solver);
assert(norm(polar.coeff - legacy_polar.coeff, 'fro') < 1e-12);

vcoul = 0.5 + rand(numel(idx_q), 1);
screened = isdf.screened_w(space, vcoul, polar);
legacy_screened = isdf_static_screened_interaction(space, vcoul, polar);
assert(norm(screened.k_mu - legacy_screened.k_mu, 'fro') < 1e-12);

target = randn(numel(idx_q), 4) + 1i * randn(numel(idx_q), 4);
kernel = isdf.screened_kernel(screened, target, vcoul);
legacy_kernel = isdf_screened_coulomb_kernel( ...
    legacy_screened, target, vcoul);
assert(norm(kernel - legacy_kernel, 'fro') < 1e-12);

coeff = randn(4, 1) + 1i * randn(4, 1);
assert(abs(isdf.screened_contract(kernel, coeff) - ...
    isdf_screened_coulomb_contract(legacy_kernel, coeff)) < 1e-12);

fprintf('ISDF package API compatibility test passed.\n');
