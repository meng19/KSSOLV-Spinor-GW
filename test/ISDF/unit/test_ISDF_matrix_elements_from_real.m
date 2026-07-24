clc;
clear;

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(31, 'twister');

fftgrid = [4, 3, 2];
nfft = prod(fftgrid);
nphi = 2;
npsi = 3;
idx_q = [1; 2; 5; 9; 13; 24];

phi = randn(nfft, nphi) + 1i * randn(nfft, nphi);
psi = randn(nfft, npsi) + 1i * randn(nfft, npsi);

direct = zeros(length(idx_q), nphi, npsi);
for iphi = 1:nphi
    for ipsi = 1:npsi
        product_grid = reshape(phi(:, iphi) .* psi(:, ipsi), fftgrid);
        product_g = fftn(product_grid) / nfft;
        direct(:, iphi, ipsi) = product_g(idx_q);
    end
end

isdf_options.rank = nphi * npsi;
isdf_options.sample_method = 'qrcp';
isdf_options.seed = 12;

actual = isdf.matrix_elements(conj(phi), psi, idx_q, fftgrid, isdf_options);

max_error = max(abs(actual(:) - direct(:)));
assert(max_error < 1e-10, ...
    'Generic ISDF matrix elements differ from direct FFT result: %.3e', max_error);

fprintf('Generic ISDF matrix element test passed. max_error = %.3e\n', max_error);
