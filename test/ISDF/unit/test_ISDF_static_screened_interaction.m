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
epsilon_vcoul = abs(randn(ng, 1)) + 0.5;
a = randn(nmu);
polar = struct();
polar.coeff = -(a * a' + eye(nmu));

screened = isdf.screened_w(vc_space, epsilon_vcoul, polar);

vmat = vc_space.zeta_g' * diag(epsilon_vcoul) * vc_space.zeta_g;
expected_eps_mu = inv(polar.coeff) - vmat;
expected_k_mu = inv(expected_eps_mu);

assert(norm(screened.smw_denominator - expected_eps_mu, 'fro') < 1e-10);
assert(norm(screened.k_mu - expected_k_mu, 'fro') < 1e-10);
assert(norm(screened.coulomb_reduced - vmat, 'fro') < 1e-10);
assert(norm(screened.polar_coeff - polar.coeff, 'fro') < 1e-10);

target_rank = 3;
target_zeta = randn(ng, target_rank) + 1i * randn(ng, target_rank);
contract_vcoul = abs(randn(ng, 1)) + 0.25;

eps_full = eye(ng) - diag(epsilon_vcoul) * vc_space.zeta_g * polar.coeff * vc_space.zeta_g';
full_screened_coul = (inv(eps_full) - eye(ng)) * diag(contract_vcoul);
actual_matrix = isdf.screened_kernel(screened, [], contract_vcoul);
assert(norm(actual_matrix - full_screened_coul, 'fro') / norm(full_screened_coul, 'fro') < 1e-10, ...
    'Reduced SMW screened Coulomb matrix does not match full inverse.');

expected_kernel = target_zeta.' * full_screened_coul * conj(target_zeta);
actual_kernel = isdf.screened_kernel(screened, target_zeta, contract_vcoul);

assert(norm(actual_kernel - expected_kernel, 'fro') / norm(expected_kernel, 'fro') < 1e-10, ...
    'Reduced SMW screened Coulomb kernel does not match full inverse.');

same_coul_full = (inv(eps_full) - eye(ng)) * diag(epsilon_vcoul);
same_coul_expected = target_zeta.' * same_coul_full * conj(target_zeta);
same_coul_actual = isdf.screened_kernel(screened, target_zeta, epsilon_vcoul);
assert(norm(same_coul_actual - same_coul_expected, 'fro') / norm(same_coul_expected, 'fro') < 1e-10, ...
    'Reduced SMW screened Coulomb kernel same-Coulomb fast path does not match full inverse.');

coeff_left = randn(target_rank, 1) + 1i * randn(target_rank, 1);
coeff_right = randn(target_rank, 1) + 1i * randn(target_rank, 1);
expected_value = coeff_left.' * expected_kernel * conj(coeff_right);
actual_value = isdf.screened_contract(actual_kernel, coeff_left, coeff_right);

assert(abs(actual_value - expected_value) / max(1, abs(expected_value)) < 1e-10, ...
    'Reduced SMW screened Coulomb contraction does not match full contraction.');

fprintf('ISDF static screened interaction test passed.\n');
