clc;
clear;

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(13, 'twister');

fftgrid = [4, 3, 2];
nfft = prod(fftgrid);
nbands = 4;
idx_q = [1; 3; 6; 11; 17; 23];

state_real = randn(nfft, 1) + 1i * randn(nfft, 1);
sum_real = randn(nfft, nbands) + 1i * randn(nfft, nbands);

direct = zeros(length(idx_q), nbands);
for nn = 1:nbands
    product_grid = reshape(conj(state_real) .* sum_real(:, nn), fftgrid);
    product_g = fftn(product_grid) / nfft;
    direct(:, nn) = product_g(idx_q);
end

isdf_options.rank = nbands;
isdf_options.sample_method = 'qrcp';
isdf_options.seed = 9;

actual = isdf_sigma_matrix_elements_from_real(state_real, sum_real, idx_q, fftgrid, isdf_options);

max_error = max(abs(actual(:) - direct(:)));
assert(max_error < 1e-10, ...
    'ISDF sigma matrix elements differ from direct FFT result: %.3e', max_error);

fprintf('ISDF sigma matrix element test passed. max_error = %.3e\n', max_error);
