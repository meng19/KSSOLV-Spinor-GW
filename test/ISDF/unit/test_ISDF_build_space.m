script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(11);
ngrid = 16;
nphi = 3;
npsi = 4;
phi = randn(ngrid, nphi) + 1i * randn(ngrid, nphi);
psi = randn(ngrid, npsi) + 1i * randn(ngrid, npsi);
idx_q = (1:ngrid).';
fftgrid = [4, 4, 1];

options = struct();
options.rank = nphi * npsi;
options.sample_method = 'qrcp';
options.seed = 0;

space = isdf.build_space(conj(phi), psi, idx_q, fftgrid, options);

assert(isfield(space, 'ind_mu'));
assert(isfield(space, 'zeta_g'));
assert(isfield(space, 'phi_mu'));
assert(isfield(space, 'psi_mu'));
assert(numel(space.ind_mu) == options.rank);
assert(size(space.zeta_g, 2) == options.rank);
assert(size(space.phi_mu, 1) == options.rank);
assert(size(space.psi_mu, 1) == options.rank);

fprintf('ISDF build space test passed.\n');
