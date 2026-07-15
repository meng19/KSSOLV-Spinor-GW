script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(13);
ng = 8;
nmu = 4;
vc_space = struct();
vc_space.zeta_g = randn(ng, nmu) + 1i * randn(ng, nmu);
vcoul = abs(randn(ng, 1)) + 0.5;
a = randn(nmu);
polar = struct();
polar.coeff = -(a * a' + eye(nmu));

screened = isdf_static_screened_interaction(vc_space, vcoul, polar);

vmat = vc_space.zeta_g' * diag(vcoul) * vc_space.zeta_g;
expected_eps_mu = inv(polar.coeff) - vmat;
expected_eps_mu_inv = inv(expected_eps_mu);
expected_w_mu = expected_eps_mu_inv;

assert(norm(screened.eps_mu - expected_eps_mu, 'fro') < 1e-10);
assert(norm(screened.eps_mu_inv - expected_eps_mu_inv, 'fro') < 1e-10);
assert(norm(screened.w_mu - expected_w_mu, 'fro') < 1e-10);

target_rank = 3;
target_zeta = randn(ng, target_rank) + 1i * randn(ng, target_rank);
right_vcoul = abs(randn(ng, 1)) + 0.25;

eps_full = eye(ng) - diag(vcoul) * vc_space.zeta_g * polar.coeff * vc_space.zeta_g';
full_screened_coul = (inv(eps_full) - eye(ng)) * diag(right_vcoul);
expected_kernel = target_zeta.' * full_screened_coul * conj(target_zeta);
actual_kernel = isdf_screened_coulomb_kernel(screened, target_zeta, right_vcoul);

assert(norm(actual_kernel - expected_kernel, 'fro') / norm(expected_kernel, 'fro') < 1e-10, ...
    'Reduced SMW screened Coulomb kernel does not match full inverse.');

coeff_left = randn(target_rank, 1) + 1i * randn(target_rank, 1);
coeff_right = randn(target_rank, 1) + 1i * randn(target_rank, 1);
expected_value = coeff_left.' * expected_kernel * conj(coeff_right);
actual_value = isdf_screened_coulomb_contract(actual_kernel, coeff_left, coeff_right);

assert(abs(actual_value - expected_value) / max(1, abs(expected_value)) < 1e-10, ...
    'Reduced SMW screened Coulomb contraction does not match full contraction.');

fprintf('ISDF static screened interaction test passed.\n');
