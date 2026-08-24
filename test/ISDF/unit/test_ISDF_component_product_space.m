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

left = cell(1, ncomponents);
right = cell(1, ncomponents);
for icomponent = 1:ncomponents
    left{icomponent} = randn(nfft, nleft) + ...
        1i * randn(nfft, nleft);
    right{icomponent} = randn(nfft, nright) + ...
        1i * randn(nfft, nright);
end

products = zeros(nfft, nleft * nright);
for iright = 1:nright
    for ileft = 1:nleft
        column = ileft + (iright - 1) * nleft;
        for icomponent = 1:ncomponents
            products(:, column) = products(:, column) + ...
                conj(left{icomponent}(:, ileft)) .* ...
                right{icomponent}(:, iright);
        end
    end
end

direct = zeros(numel(idx_q), size(products, 2));
for ipair = 1:size(products, 2)
    product_grid = reshape(products(:, ipair), fftgrid);
    product_g = fftn(product_grid) / nfft;
    direct(:, ipair) = product_g(idx_q);
end

options.rank = nleft * nright;
options.sample_method = 'qrcp';
options.seed = 0;
space = isdf.build_space(left, right, idx_q, fftgrid, options);
actual = space.zeta_g * space.product_mu;
max_error = max(abs(actual(:) - direct(:)));
assert(max_error < 1e-10, ...
    'Component product-space ISDF differs from direct FFT: %.3e', ...
    max_error);
selected_products = products(space.ind_mu, :);
assert(max(abs(space.product_mu(:) - selected_products(:))) < 1e-12);

ev_occ = [-0.9; -0.2];
ev_unocc = [0.4; 0.8; 1.5];
direct_polar = isdf.polarizability(space, ev_occ, ev_unocc, ...
    struct('method', 'direct'));
cauchy_polar = isdf.polarizability(space, ev_occ, ev_unocc, ...
    struct('method', 'cauchy', 'froErr', 1e-10, 'MaxIter', 10));
relative_error = norm(direct_polar.coeff - cauchy_polar.coeff, 'fro') / ...
    max(1, norm(direct_polar.coeff, 'fro'));
assert(relative_error < 1e-8, ...
    'Component Cauchy coefficient mismatch: %.3e', relative_error);

single_products = zeros(nfft, nleft * nright);
for iright = 1:nright
    for ileft = 1:nleft
        column = ileft + (iright - 1) * nleft;
        single_products(:, column) = conj(left{1}(:, ileft)) .* ...
            right{1}(:, iright);
    end
end
single_space = isdf.build_space(left{1}, right{1}, ...
    idx_q, fftgrid, options);
assert(isequal(single_space.left_mu_components{1}, ...
    left{1}(single_space.ind_mu, :)));
assert(isequal(single_space.right_mu_components{1}, ...
    right{1}(single_space.ind_mu, :)));
single_direct = zeros(numel(idx_q), size(single_products, 2));
for ipair = 1:size(single_products, 2)
    product_g = fftn(reshape(single_products(:, ipair), fftgrid)) / nfft;
    single_direct(:, ipair) = product_g(idx_q);
end
single_actual = single_space.zeta_g * single_space.product_mu;
single_error = max(abs(single_actual(:) - single_direct(:)));
assert(single_error < 1e-10, ...
    'Single-component product-space mismatch: %.3e', single_error);

randomized_options = options;
randomized_options.sample_method = 'qrcp_randomized';
randomized_options.random_oversampling = 1;
randomized_space = isdf.build_space(left, right, ...
    idx_q, fftgrid, randomized_options);
randomized_actual = randomized_space.zeta_g * randomized_space.product_mu;
randomized_error = max(abs(randomized_actual(:) - direct(:)));
assert(randomized_error < 1e-10, ...
    'Randomized component ISDF differs from direct FFT: %.3e', ...
    randomized_error);
assert(~isfield(randomized_space, 'products'));

kmeans_options = options;
kmeans_options.sample_method = 'kmeans';
kmeans_options.kmeans_max_iter = 4;
kmeans_options.kmeans_replicates = 1;
kmeans_space = isdf.build_space(left, right, ...
    idx_q, fftgrid, kmeans_options);
kmeans_actual = kmeans_space.zeta_g * kmeans_space.product_mu;
kmeans_error = max(abs(kmeans_actual(:) - direct(:)));
assert(kmeans_error < 1e-10, ...
    'K-means component ISDF differs from direct FFT: %.3e', kmeans_error);
assert(~isfield(kmeans_space, 'products'));

explicit_add_options = kmeans_options;
explicit_add_options.weight = 'add';
explicit_add_space = isdf.build_space(left, right, ...
    idx_q, fftgrid, explicit_add_options);
assert(numel(explicit_add_space.ind_mu) == explicit_add_options.rank);

explicit_power_options = kmeans_options;
explicit_power_options.weight = 'power';
explicit_power_options.power = 1.5;
explicit_power_space = isdf.build_space(left, right, ...
    idx_q, fftgrid, explicit_power_options);
assert(numel(explicit_power_space.ind_mu) == explicit_power_options.rank);

bad_weight_options = kmeans_options;
bad_weight_options.weight = 'bad';
try
    isdf.build_space(left, right, idx_q, fftgrid, bad_weight_options);
    error('ISDFTest:MissingUnknownWeightError', ...
        'Component K-means accepted an invalid weight option.');
catch err
    assert(strcmp(err.identifier, 'ISDF:UnknownWeight'), ...
        'Unexpected invalid component weight error: %s', err.identifier);
end

fprintf(['Component package ISDF tests passed. qrcp = %.3e, ' ...
    'single = %.3e, randomized = %.3e, kmeans = %.3e, cauchy = %.3e\n'], ...
    max_error, single_error, randomized_error, kmeans_error, relative_error);

% A single-component QRCP space must retain the legacy Gram solve.  Its
% condition estimate is for P_mu*P_mu', rather than for rectangular P_mu.
rank_left = randn(12, 2) + 1i * randn(12, 2);
rank_right = randn(12, 2) + 1i * randn(12, 2);
rank_options = struct('rank', 3, 'sample_method', 'qrcp', ...
    'seed', 4, 'rcond_tol', 1, ...
    'warn_ill_conditioned', false, 'ill_conditioned_solver', 'pinv');
rank_space = isdf.build_space(rank_left, rank_right, (1:12).', ...
    [3, 4, 1], rank_options);

rank_products = zeros(12, 4);
for iright = 1:2
    for ileft = 1:2
        column = ileft + (iright - 1) * 2;
        rank_products(:, column) = conj(rank_left(:, ileft)) .* ...
            rank_right(:, iright);
    end
end
rank_products_mu = rank_products(rank_space.ind_mu, :);
rank_c1 = rank_products * rank_products_mu';
rank_c2 = rank_products_mu * rank_products_mu';
rank_zeta_real = rank_c1 * pinv(rank_c2);
rank_zeta_reference = zeros(12, rank_options.rank);
for imu = 1:rank_options.rank
    rank_fft = fftn(reshape(rank_zeta_real(:, imu), [3, 4, 1])) / 12;
    rank_zeta_reference(:, imu) = rank_fft(:);
end
assert(abs(rank_space.solve_info.rcond - rcond(rank_c2)) < 1e-15, ...
    'Single-component QRCP solve_info must describe the Gram denominator.');
assert(rank_space.solve_info.used_pinv);
assert(max(abs(rank_space.zeta_g(:) - rank_zeta_reference(:))) < 1e-10, ...
    'Single-component QRCP must use the locked Gram pseudoinverse path.');

deficient_left = rank_left;
deficient_left(:, 2) = deficient_left(:, 1);
deficient_options = rmfield(rank_options, 'ill_conditioned_solver');
deficient_options.rcond_tol = 1e-12;
warning_state = warning('query', 'MATLAB:rankDeficientMatrix');
warning_cleanup = onCleanup(@() warning( ...
    warning_state.state, 'MATLAB:rankDeficientMatrix'));
warning('on', 'MATLAB:rankDeficientMatrix');
lastwarn('');
isdf.build_space(deficient_left, rank_right, (1:12).', ...
    [3, 4, 1], deficient_options);
[~, warning_id] = lastwarn;
assert(~strcmp(warning_id, 'MATLAB:rankDeficientMatrix'), ...
    'Single-component Gram solve leaked a rank-deficient warning.');

% Structural memory gates: component K-means weights stay O(ngrid), and
% explicit QRCP products are returned to build_space instead of rebuilt.
sample_source = fileread(fullfile(repo_root, ...
    'src', 'GW', 'ISDF', '+isdf', 'private', 'sample_points.m'));
weight_source = fileread(fullfile(repo_root, ...
    'src', 'GW', 'ISDF', '+isdf', 'private', 'component_weight.m'));
build_source = fileread(fullfile(repo_root, ...
    'src', 'GW', 'ISDF', '+isdf', 'build_space.m'));
assert(contains(weight_source, 'function weight = component_weight'));
assert(contains(weight_source, 'switch lower(options.weight)'), ...
    'Component K-means weight helper must honor options.weight.');
assert(~contains(weight_source, 'component_products('), ...
    'Component K-means weight helper must not form the full product matrix.');
assert(contains(build_source, ...
    '[ind_mu, products] = sample_points(left, right, options);'), ...
    'build_space must reuse explicit products returned by sample_points.');
explicit_builds = regexp(sample_source, ...
    'component_products\(left, right, \[\], \[\]\)', 'match');
assert(numel(explicit_builds) == 1, ...
    'The explicit component product matrix must be formed only once.');
