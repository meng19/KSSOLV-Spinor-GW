script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

if exist('gpuDevice', 'file') ~= 2 || gpuDeviceCount < 1
    fprintf('ISDF GPU support test skipped: no GPU device available.\n');
    return;
end

gpuDevice(1);
rng(17);

fftgrid = [4, 4, 1];
ngrid = prod(fftgrid);
nleft = 2;
nright = 3;
nmu = 3;
idx_q = gpuArray((1:ngrid).');

left = cell(1, 2);
right = cell(1, 2);
for icomponent = 1:2
    left{icomponent} = gpuArray( ...
        randn(ngrid, nleft) + 1i * randn(ngrid, nleft));
    right{icomponent} = gpuArray( ...
        randn(ngrid, nright) + 1i * randn(ngrid, nright));
end

options = struct();
options.rank = nmu;
options.sample_method = 'qrcp_randomized';
options.seed = 17;
options.fftgrid = fftgrid;
options.rcond_tol = 1e-12;
options.warn_ill_conditioned = false;

space = isdf.build_space(left, right, idx_q, fftgrid, options);
assert(isa(space.zeta_g, 'gpuArray'));
assert(isa(space.product_mu, 'gpuArray'));

solver = struct();
solver.method = 'direct';
solver.freq = 0;
ev_occ = [-2.0; -1.0];
ev_unocc = [0.5; 1.5; 2.5];
polar = isdf.polarizability(space, ev_occ, ev_unocc, solver);
assert(isa(polar.coeff, 'gpuArray'));

epsilon_vcoul = gpuArray(abs(randn(ngrid, 1)) + 0.5);
screened = isdf.screened_w(space, epsilon_vcoul, polar);
assert(isa(screened.k_mu, 'gpuArray'));

kernel = isdf.screened_kernel(screened, space.zeta_g, epsilon_vcoul);
assert(isa(kernel, 'gpuArray'));

value = isdf.screened_contract(kernel(:, :, 1), space.product_mu(:, 1));
assert(isa(value, 'gpuArray'));

gme = isdf.matrix_elements(left, right, idx_q, fftgrid, options);
assert(isa(gme, 'gpuArray'));
assert(all(size(gme) == [ngrid, nleft, nright]));

reset(gpuDevice);
fprintf('ISDF GPU support smoke test passed.\n');
