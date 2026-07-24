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
subset_products = isdf_component_products(left_components, right_components, idx_q);
subset_products_ref = products(idx_q, :);
assert(max(abs(subset_products(:) - subset_products_ref(:))) < 1e-12, ...
    'Component products with grid_indices should match selected full-product rows.');

single_subset_products = isdf_component_products({left_components{1}}, ...
    {right_components{1}}, idx_q);
single_subset_ref = isdf_prod_states(conj(left_components{1}(idx_q, :)), ...
    right_components{1}(idx_q, :));
assert(max(abs(single_subset_products(:) - single_subset_ref(:))) < 1e-12, ...
    'Single-component indexed products should match scalar product states.');

component_weight = isdf_component_product_weight(left_components, right_components);
component_weight_ref = sum(abs(products).^2, 2);
assert(max(abs(component_weight(:) - component_weight_ref(:))) < 1e-10, ...
    'Component product weight should match explicit product row norms.');

projection = randn(nleft * nright, 4) + 1i * randn(nleft * nright, 4);
projected_products = isdf_component_products_projected( ...
    left_components, right_components, projection);
projected_products_ref = products * projection;
assert(max(abs(projected_products(:) - projected_products_ref(:))) < 1e-10, ...
    'Projected component products should match explicit product projection.');

direct = zeros(length(idx_q), nleft * nright);
for ipair = 1:size(products, 2)
    product_grid = reshape(products(:, ipair), fftgrid);
    product_g = fftn(product_grid) / nfft;
    direct(:, ipair) = product_g(idx_q);
end

ind_mu_ref = [2; 7; 11; 18; 21; 24];
[component_c1, component_c2] = isdf_product_gram( ...
    left_components, right_components, ind_mu_ref);
component_products_mu = products(ind_mu_ref, :);
component_c1_ref = products * component_products_mu';
component_c2_ref = component_products_mu * component_products_mu';
assert(max(abs(component_c1(:) - component_c1_ref(:))) < 1e-12, ...
    'Component product Gram c1 should match explicit product matrix.');
assert(max(abs(component_c2(:) - component_c2_ref(:))) < 1e-12, ...
    'Component product Gram c2 should match explicit product matrix.');

scalar_products = isdf_prod_states(conj(left_components{1}), right_components{1});
[scalar_c1, scalar_c2] = isdf_product_gram( ...
    {left_components{1}}, {right_components{1}}, ind_mu_ref);
scalar_products_mu = scalar_products(ind_mu_ref, :);
scalar_c1_ref = scalar_products * scalar_products_mu';
scalar_c2_ref = scalar_products_mu * scalar_products_mu';
assert(max(abs(scalar_c1(:) - scalar_c1_ref(:))) < 1e-12, ...
    'Single-component product Gram c1 should match explicit product matrix.');
assert(max(abs(scalar_c2(:) - scalar_c2_ref(:))) < 1e-12, ...
    'Single-component product Gram c2 should match explicit product matrix.');

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

scalar_product_mu = isdf_prod_states(scalar_space.phi_mu, scalar_space.psi_mu);
assert(max(abs(single_component_space.product_mu(:) - scalar_product_mu(:))) < 1e-12, ...
    'Single-component product_mu should match scalar product_mu.');

single_component_reconstruction = single_component_space.zeta_g * single_component_space.product_mu;
scalar_reconstruction = scalar_space.zeta_g * scalar_product_mu;
single_component_error = max(abs(single_component_reconstruction(:) - scalar_reconstruction(:)));
assert(single_component_error < 1e-10, ...
    'Single-component space should reuse scalar ISDF representation: %.3e', ...
    single_component_error);

fprintf('Single-component fast-path test passed. max_error = %.3e\n', ...
    single_component_error);

matrix_free_options = isdf_options;
matrix_free_options.sample_method = 'qrcp_randomized';
matrix_free_options.random_oversampling = 1;
matrix_free_space = isdf_build_space(left_components, right_components, ...
    idx_q, fftgrid, matrix_free_options);
matrix_free_reconstruction = matrix_free_space.zeta_g * matrix_free_space.product_mu;
matrix_free_error = max(abs(matrix_free_reconstruction(:) - direct(:)));
assert(matrix_free_error < 1e-10, ...
    'Matrix-free randomized component ISDF differs from direct FFT result: %.3e', ...
    matrix_free_error);
assert(~isfield(matrix_free_space, 'products'), ...
    'Matrix-free randomized component ISDF should not store explicit products.');

kmeans_options = isdf_options;
kmeans_options.sample_method = 'kmeans';
kmeans_options.kmeans_max_iter = 4;
kmeans_options.kmeans_replicates = 1;
kmeans_space = isdf_build_space(left_components, right_components, ...
    idx_q, fftgrid, kmeans_options);
kmeans_reconstruction = kmeans_space.zeta_g * kmeans_space.product_mu;
kmeans_error = max(abs(kmeans_reconstruction(:) - direct(:)));
assert(kmeans_error < 1e-10, ...
    'Matrix-free kmeans component ISDF differs from direct FFT result: %.3e', ...
    kmeans_error);
assert(~isfield(kmeans_space, 'products'), ...
    'Matrix-free kmeans component ISDF should not store explicit products.');

fprintf('Matrix-free component product-space tests passed. randomized_error = %.3e, kmeans_error = %.3e\n', ...
    matrix_free_error, kmeans_error);
