% TEST_ISDF_LEFT_RIGHT_SAMPLING_SWAP Test per-space swapped point selection.
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(genpath(fullfile(root, 'src')));

rng(20260921);
fftgrid = [4, 4, 2];
ngrid = prod(fftgrid);
left = randn(ngrid, 5) + 1i * randn(ngrid, 5);
right = randn(ngrid, 7) + 1i * randn(ngrid, 7);
options = struct('rank', 8, 'sample_method', 'qrcp_randomized', ...
    'random_oversampling', 2, 'seed', 11, ...
    'warn_rank_selection', false, 'interpolation_solver', 'direct');
reference = isdf.build_space(left, right, (1:ngrid).', fftgrid, options);
swapped_options = options;
swapped_options.swap_left_right = true;
swapped = isdf.build_space(left, right, (1:ngrid).', fftgrid, ...
    swapped_options);
assert(swapped.sampling_swapped_left_right);
assert(isequal(size(reference.product_mu), size(swapped.product_mu)), ...
    'Sampling swap must retain the original physical product dimensions.');
assert(isequal(size(reference.zeta_g), size(swapped.zeta_g)), ...
    'Sampling swap must retain the original output dimensions.');
assert(~isequal(reference.ind_mu, swapped.ind_mu), ...
    'Randomized left-right sampling control unexpectedly chose identical pivots.');
fprintf('test_ISDF_left_right_sampling_swap passed\n');
