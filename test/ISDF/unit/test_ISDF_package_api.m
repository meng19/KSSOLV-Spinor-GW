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

ev_occ = [-0.8; -0.2];
ev_unocc = [0.3; 0.9; 1.4];
solver = struct('method', 'direct');
polar = isdf.polarizability(space, ev_occ, ev_unocc, solver);
polar_reference = zeros(space.rank);
for iv = 1:numel(ev_occ)
    for ic = 1:numel(ev_unocc)
        coefficient = space.product_mu(:, iv + ...
            (ic - 1) * numel(ev_occ));
        polar_reference = polar_reference + coefficient * coefficient' / ...
            (ev_occ(iv) - ev_unocc(ic));
    end
end
assert(norm(polar.coeff - polar_reference, 'fro') < 1e-12);

vcoul = 0.5 + rand(numel(idx_q), 1);
screened = isdf.screened_w(space, vcoul, polar);
reduced_coulomb = space.zeta_g' * (vcoul .* space.zeta_g);
screened_reference = (eye(size(polar.coeff)) - ...
    polar.coeff * reduced_coulomb) \ polar.coeff;
assert(norm(screened.k_mu - screened_reference, 'fro') < 1e-12);

target = randn(numel(idx_q), 4) + 1i * randn(numel(idx_q), 4);
kernel = isdf.screened_kernel(screened, target, vcoul);
left_projector = target.' * (vcoul .* space.zeta_g);
kernel_reference = left_projector * screened.k_mu * left_projector';
assert(norm(kernel - kernel_reference, 'fro') < 1e-12);

coeff = randn(4, 1) + 1i * randn(4, 1);
contract_reference = coeff.' * kernel * conj(coeff);
assert(abs(isdf.screened_contract(kernel, coeff) - ...
    contract_reference) < 1e-12);

fprintf('ISDF package API compatibility test passed.\n');

prefix = 'isdf_';
flat_names = { ...
    [prefix 'build_space'], ...
    [prefix 'reduced_polarizability'], ...
    [prefix 'static_screened_interaction'], ...
    [prefix 'screened_coulomb_kernel'], ...
    [prefix 'screened_coulomb_contract'], ...
    [prefix 'epsilon_batch'], ...
    [prefix 'sigma_batch'], ...
    [prefix 'wavefunction_real_component'], ...
    [prefix 'comega_cstar']};
for i = 1:numel(flat_names)
    assert(exist(flat_names{i}, 'file') == 0, ...
        'Superseded flat ISDF API remains: %s', flat_names{i});
end

assert(~isempty(which('isdf.build_space')));
assert(~isempty(which('isdf.matrix_elements')));
