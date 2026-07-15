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
ngrid = prod(fftgrid);
nphi = 3;
npsi = 4;
rank_mu = 6;

phi = randn(ngrid, nphi) + 1i * randn(ngrid, nphi);
psi = randn(ngrid, npsi) + 1i * randn(ngrid, npsi);
idx_q = (1:ngrid).';

methods = {'qrcp_randomized', 'kmeans', 'default'};
for imethod = 1:numel(methods)
    options = struct();
    options.rank = rank_mu;
    options.sample_method = methods{imethod};
    options.seed = 5;
    options.weight = 'add';
    options.init = 'kmeans++';
    options.kmeans_max_iter = 20;
    options.fftgrid = fftgrid;

    ind_mu = isdf_indices(phi, psi, options);
    assert(numel(ind_mu) == rank_mu, ...
        '%s returned %d points instead of %d.', methods{imethod}, numel(ind_mu), rank_mu);
    assert(numel(unique(ind_mu)) == rank_mu, ...
        '%s returned duplicate interpolation points.', methods{imethod});
    assert(all(ind_mu >= 1) && all(ind_mu <= ngrid), ...
        '%s returned out-of-range interpolation points.', methods{imethod});
end

options = struct();
options.rank = nphi * npsi;
options.sample_method = 'qrcp_randomized';
options.seed = 9;
options.fftgrid = fftgrid;

actual = isdf_matrix_elements_from_real(phi, psi, idx_q, fftgrid, options);
direct = zeros(length(idx_q), nphi, npsi);
for iphi = 1:nphi
    for ipsi = 1:npsi
        product_grid = reshape(phi(:, iphi) .* psi(:, ipsi), fftgrid);
        product_g = fftn(product_grid) / ngrid;
        direct(:, iphi, ipsi) = product_g(idx_q);
    end
end

max_error = max(abs(actual(:) - direct(:)));
assert(max_error < 1e-10, ...
    'Randomized QRCP full-rank matrix elements differ from direct FFT: %.3e', max_error);

fprintf('ISDF interpolation point method test passed. max_error = %.3e\n', max_error);
