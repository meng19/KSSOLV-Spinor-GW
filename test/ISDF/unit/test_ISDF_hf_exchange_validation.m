% TEST_ISDF_HF_EXCHANGE_VALIDATION Validate the direct/ISDF EX comparison.
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(genpath(fullfile(root, 'src')));

rng(17);
fftgrid = [2, 2, 1];
ngrid = prod(fftgrid);
left = randn(ngrid, 2) + 1i * randn(ngrid, 2);
right = randn(ngrid, 2) + 1i * randn(ngrid, 2);
options = struct('rank', 4, 'sample_method', 'qrcp', ...
    'warn_rank_selection', false);
space = isdf.build_space(left, right, (1:ngrid).', fftgrid, options);
weights = [1; 0.5; 0.25; 0];
report = isdf.validate_hf_exchange(space, left, right, (1:ngrid).', ...
    fftgrid, (1:ngrid).', weights, struct());
assert(report.pairs_checked == 3);
assert(report.relative_error < 1e-11, ...
    'Full-rank ISDF must reproduce the direct HF exchange energy.');

limited = isdf.validate_hf_exchange(space, left, right, (1:ngrid).', ...
    fftgrid, (1:ngrid).', weights, ...
    struct('validate_hf_exchange_max_pairs', 2));
assert(limited.pairs_checked == 2);
fprintf('test_ISDF_hf_exchange_validation passed\n');
