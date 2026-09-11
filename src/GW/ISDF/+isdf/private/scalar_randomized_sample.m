function ind_mu = scalar_randomized_sample(left, right, options)
%SCALAR_RANDOMIZED_SAMPLE Select scalar ISDF points via randomized QRCP.

nleft = size(left, 2);
nright = size(right, 2);
rank_mu = options.rank;
sample_rank = max(rank_mu, ...
    ceil(options.random_oversampling * rank_mu));
left_rank = min(max(1, ceil(sqrt((nleft / nright) * sample_rank))), ...
    nleft);
right_rank = min(max(1, ceil(sqrt((nright / nleft) * sample_rank))), ...
    nright);

% This sketch is used only to select interpolation points. Keep the final
% product values and interpolation solve in the wavefunction's native
% precision; sample_precision may reduce only this temporary QRCP.
local_progress(options, 0.04, 'sampling randomized sketch');
left_sample = sample_cast(left, options);
right_sample = sample_cast(right, options);
left_projection = randn_like(nleft, left_rank, left_sample);
right_projection = randn_like(nright, right_rank, right_sample);
if ~isreal(left) || ~isreal(right)
    left_projection = left_projection + ...
        1i * randn_like(nleft, left_rank, left);
    right_projection = right_projection + ...
        1i * randn_like(nright, right_rank, right);
end
local_progress(options, 0.08, 'sampling left projection');
compressed_left = conj(left_sample) * left_projection;
local_progress(options, 0.14, 'sampling right projection');
compressed_right = right_sample * right_projection;
local_progress(options, 0.20, 'sampling sketch products');
products = pair_products(compressed_left, compressed_right, ...
    @(current, total) local_product_progress(options, current, total));
local_progress(options, 0.30, 'sampling QRCP');
ind_mu = qrcp_sample(products, rank_mu);
end

function local_product_progress(options, current, total)
fraction = 0.20 + 0.08 * current / total;
local_progress(options, fraction, sprintf( ...
    'sampling sketch products %d/%d', current, total));
end

function local_progress(options, fraction, stage)
if isfield(options, 'progress') && isa(options.progress, 'function_handle')
    options.progress(fraction, stage);
end
end
