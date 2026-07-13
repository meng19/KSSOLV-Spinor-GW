clc;
clear;

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(11, 'twister');

fftgrid = [4, 3, 2];
nfft = prod(fftgrid);
nv = 2;
nc = 3;
idx_q = [1; 2; 5; 9; 13; 24];

valence_real = randn(nfft, nv) + 1i * randn(nfft, nv);
conduction_real = randn(nfft, nc) + 1i * randn(nfft, nc);

phi = conj(valence_real);
psi = conduction_real;

direct = zeros(length(idx_q), nv, nc);
for iv = 1:nv
    for ic = 1:nc
        product_grid = reshape(phi(:, iv) .* psi(:, ic), fftgrid);
        product_g = fftn(product_grid) / nfft;
        direct(:, iv, ic) = product_g(idx_q);
    end
end

isdf_options.rank = nv * nc;
isdf_options.sample_method = 'qrcp';
isdf_options.seed = 7;

actual = isdf_epsilon_matrix_elements_from_real(phi, psi, idx_q, fftgrid, isdf_options);

max_error = max(abs(actual(:) - direct(:)));
assert(max_error < 1e-10, ...
    'ISDF epsilon matrix elements differ from direct FFT result: %.3e', max_error);

fprintf('ISDF epsilon matrix element test passed. max_error = %.3e\n', max_error);
