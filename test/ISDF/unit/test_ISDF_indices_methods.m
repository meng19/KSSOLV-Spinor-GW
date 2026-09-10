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

typed_options = struct('rank', 99, 'rank_vc', 5, ...
    'rank_ratio_nn', 2);
vc_options = isdf.options_for_type(typed_options, 'vc');
assert(vc_options.rank == 5, 'rank_vc must override the generic rank.');
nn_options = isdf.options_for_type(typed_options, 'nn');
assert(~isfield(nn_options, 'rank') && nn_options.rank_ratio == 2, ...
    'rank_ratio_nn must request automatic NN rank selection.');

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

    space = isdf.build_space(conj(phi), psi, idx_q, fftgrid, options);
    ind_mu = space.ind_mu;
    assert(numel(ind_mu) == rank_mu, ...
        '%s returned %d points instead of %d.', methods{imethod}, numel(ind_mu), rank_mu);
    assert(numel(unique(ind_mu)) == rank_mu, ...
        '%s returned duplicate interpolation points.', methods{imethod});
    assert(all(ind_mu >= 1) && all(ind_mu <= ngrid), ...
        '%s returned out-of-range interpolation points.', methods{imethod});
end

adaptive_options = struct();
adaptive_options.rank = 2;
adaptive_options.sample_method = 'qrcp';
adaptive_options.adaptive_rank_enable = true;
adaptive_options.adaptive_rank_tol = 1e-12;
adaptive_options.adaptive_rank_step = 2;
adaptive_options.adaptive_rank_max = nphi * npsi;
adaptive_options.warn_rank_selection = false;
adaptive_space = isdf.build_space(conj(phi), psi, idx_q, fftgrid, ...
    adaptive_options);
assert(adaptive_space.adaptive_info.enabled, ...
    'Adaptive QRCP metadata was not recorded.');
assert(adaptive_space.adaptive_info.reached_tolerance, ...
    'Adaptive QRCP did not reach the requested residual tolerance.');
assert(adaptive_space.rank >= adaptive_options.rank && ...
    adaptive_space.rank <= adaptive_options.adaptive_rank_max, ...
    'Adaptive QRCP returned a rank outside its requested bounds.');
assert(strcmp(adaptive_space.adaptive_info.residual_kind, 'exact_residual'));

adaptive_randomized_options = adaptive_options;
adaptive_randomized_options.sample_method = 'qrcp_randomized';
adaptive_randomized_options.adaptive_validation_rank = nphi * npsi;
adaptive_randomized_space = isdf.build_space(conj(phi), psi, idx_q, ...
    fftgrid, adaptive_randomized_options);
assert(adaptive_randomized_space.adaptive_info.enabled && ...
    adaptive_randomized_space.adaptive_info.reached_tolerance, ...
    'Adaptive randomized QRCP did not reach the requested validation tolerance.');
assert(strcmp(adaptive_randomized_space.adaptive_info.residual_kind, ...
    'randomized_validation'));

options = struct();
options.rank = nphi * npsi;
options.sample_method = 'qrcp_randomized';
options.seed = 9;
options.fftgrid = fftgrid;

actual = isdf.matrix_elements(conj(phi), psi, idx_q, fftgrid, options);
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
