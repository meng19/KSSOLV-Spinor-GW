clc;
clear;

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(23, 'twister');

fftgrid = [4, 3, 2];
nfft = prod(fftgrid);
nleft = 2;
nright = 3;
ncomponents = 2;
idx_q = [1; 2; 5; 9; 13; 24];

left_components = cell(1, ncomponents);
right_components = cell(1, ncomponents);
for icomponent = 1:ncomponents
    left_components{icomponent} = randn(nfft, nleft) + 1i * randn(nfft, nleft);
    right_components{icomponent} = randn(nfft, nright) + 1i * randn(nfft, nright);
end

products = isdf_component_products(left_components, right_components);
direct = zeros(length(idx_q), nleft * nright);
for ipair = 1:size(products, 2)
    product_grid = reshape(products(:, ipair), fftgrid);
    product_g = fftn(product_grid) / nfft;
    direct(:, ipair) = product_g(idx_q);
end

isdf_options.rank = nleft * nright;
isdf_options.sample_method = 'qrcp';
isdf_options.seed = 0;
space = isdf_build_space(left_components, right_components, idx_q, fftgrid, isdf_options);
actual = space.zeta_g * space.product_mu;

max_error = max(abs(actual(:) - direct(:)));
assert(max_error < 1e-10, ...
    'Component product-space ISDF differs from direct FFT result: %.3e', max_error);

ev_occ = [-0.9; -0.2];
ev_unocc = [0.4; 0.8; 1.5];
direct_options.method = 'direct';
direct_coeff = isdf_comega_cstar(space.left_mu_components, ...
    space.right_mu_components, ev_occ, ev_unocc, direct_options);

cauchy_options.method = 'cauchy';
cauchy_options.froErr = 1e-10;
cauchy_options.MaxIter = 10;
[cauchy_coeff, cauchy_info] = isdf_comega_cstar(space.left_mu_components, ...
    space.right_mu_components, ev_occ, ev_unocc, cauchy_options);

relative_error = norm(direct_coeff - cauchy_coeff, 'fro') / ...
    max(1, norm(direct_coeff, 'fro'));
assert(relative_error < 1e-8, ...
    'Component Cauchy coefficient mismatch: %.3e', relative_error);

fprintf('Component product-space ISDF test passed. max_error = %.3e, cauchy_error = %.3e, iterations = %d\n', ...
    max_error, relative_error, cauchy_info.iterations);

scalar_phi = conj(left_components{1});
scalar_psi = right_components{1};
scalar_space = isdf_build_space(scalar_phi, scalar_psi, idx_q, fftgrid, isdf_options);
single_component_space = isdf_build_space({left_components{1}}, {right_components{1}}, ...
    idx_q, fftgrid, isdf_options);

assert(isfield(single_component_space, 'phi_mu'));
assert(isfield(single_component_space, 'psi_mu'));
assert(isequal(single_component_space.left_mu_components{1}, left_components{1}(scalar_space.ind_mu, :)));
assert(isequal(single_component_space.right_mu_components{1}, scalar_space.psi_mu));
assert(isequal(conj(single_component_space.left_mu_components{1}), scalar_space.phi_mu));

scalar_products = isdf_prod_states(scalar_space.phi_mu, scalar_space.psi_mu);
assert(max(abs(single_component_space.product_mu(:) - scalar_products(:))) < 1e-12, ...
    'Single-component product_mu should match scalar product_mu.');

single_component_reconstruction = single_component_space.zeta_g * single_component_space.product_mu;
scalar_reconstruction = scalar_space.zeta_g * scalar_products;
single_component_error = max(abs(single_component_reconstruction(:) - scalar_reconstruction(:)));
assert(single_component_error < 1e-10, ...
    'Single-component space should reuse scalar ISDF representation: %.3e', ...
    single_component_error);

fprintf('Single-component fast-path test passed. max_error = %.3e\n', ...
    single_component_error);
