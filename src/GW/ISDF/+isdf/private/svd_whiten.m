function [zeta_g, coeff, transform, info] = svd_whiten( ...
        c1, c2, product_mu, idx_q, fftgrid, options)
%SVD_WHITEN Build a truncated, spectrally whitened ISDF basis.
%   With C = U*Lambda*U', use R = U*Lambda^(-ratio) and
%   L = Lambda^(ratio-1)*U', so R*L is the truncated C^{-1}.

if ~isfield(options, 'svd_cutoff') || isempty(options.svd_cutoff)
    options.svd_cutoff = 0;
end
if ~isfield(options, 'svd_ratio') || isempty(options.svd_ratio)
    options.svd_ratio = 0.5;
end
c2_host = gather_if_gpu(c2);
c2_host = (c2_host + c2_host') / 2;
[vectors, values] = eig(c2_host, 'vector');
[values, order] = sort(real(values), 'descend');
vectors = vectors(:, order);
if isempty(values) || values(1) <= 0
    error('ISDF:SVDEmptyBasis', 'The ISDF product Gram matrix has no positive mode.');
end
keep = values > values(1) * options.svd_cutoff;
vectors = vectors(:, keep);
values = values(keep);
if isempty(values)
    error('ISDF:SVDCutoff', ...
        'svd_cutoff=%g removed every ISDF interpolation mode.', options.svd_cutoff);
end

right_factor = vectors .* reshape(values .^ (-options.svd_ratio), 1, []);
transform = reshape(values .^ (options.svd_ratio - 1), [], 1) .* vectors';
if isa(c1, 'gpuArray')
    right_factor = gpuArray(right_factor);
    transform = gpuArray(transform);
end
zeta_real = c1 * right_factor;
coeff = transform * product_mu;
[zeta_g, solve_info] = zeta_to_g(zeta_real, [], idx_q, fftgrid, options);
info = struct('method', 'svd_whiten', 'rank_before', size(c2, 1), ...
    'rank_after', numel(values), 'cutoff', options.svd_cutoff, ...
    'ratio', options.svd_ratio, 'eigenvalues', values, ...
    'rcond', values(end) / values(1), 'used_pinv', true, ...
    'ill_conditioned', false, 'zeta_info', solve_info);
end
