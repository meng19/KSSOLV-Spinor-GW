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
nspinor = 2;
idx_q = [1; 2; 5; 9; 13; 24];

left_components = cell(1, nspinor);
right_components = cell(1, nspinor);
for ispinor = 1:nspinor
    left_components{ispinor} = randn(nfft, nleft) + 1i * randn(nfft, nleft);
    right_components{ispinor} = randn(nfft, nright) + 1i * randn(nfft, nright);
end

products = isdf_spinor_products_from_real(left_components, right_components);
direct = zeros(length(idx_q), nleft * nright);
for ipair = 1:size(products, 2)
    product_grid = reshape(products(:, ipair), fftgrid);
    product_g = fftn(product_grid) / nfft;
    direct(:, ipair) = product_g(idx_q);
end

isdf_options.rank = nleft * nright;
isdf_options.sample_method = 'qrcp';
isdf_options.seed = 0;
space = isdf_spinor_build_space(left_components, right_components, idx_q, fftgrid, isdf_options);
actual = space.zeta_g * space.product_mu;

max_error = max(abs(actual(:) - direct(:)));
assert(max_error < 1e-10, ...
    'Spinor product-space ISDF differs from direct FFT result: %.3e', max_error);

ev_occ = [-0.9; -0.2];
ev_unocc = [0.4; 0.8; 1.5];
direct_options.method = 'direct';
direct_coeff = isdf_spinor_comega_cstar(space.left_mu_components, ...
    space.right_mu_components, ev_occ, ev_unocc, direct_options);

cauchy_options.method = 'cauchy';
cauchy_options.froErr = 1e-10;
cauchy_options.MaxIter = 10;
[cauchy_coeff, cauchy_info] = isdf_spinor_comega_cstar(space.left_mu_components, ...
    space.right_mu_components, ev_occ, ev_unocc, cauchy_options);

relative_error = norm(direct_coeff - cauchy_coeff, 'fro') / ...
    max(1, norm(direct_coeff, 'fro'));
assert(relative_error < 1e-8, ...
    'Spinor Cauchy coefficient mismatch: %.3e', relative_error);

fprintf('Spinor product-space ISDF test passed. max_error = %.3e, cauchy_error = %.3e, iterations = %d\n', ...
    max_error, relative_error, cauchy_info.iterations);
