function ind_mu = randomized_sample(left, right, options)
%RANDOMIZED_SAMPLE Select ISDF points from a separable random product sketch.
%
% The same left/right random projections are used for every spinor
% component.  Summing their projected products represents
% sum_s conj(left_s).*right_s without materializing a dense
% (nleft*nright)-by-nprojection projection matrix.

if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
if numel(left) ~= numel(right)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end

nleft = size(left{1}, 2);
nright = size(right{1}, 2);
rank_mu = options.rank;
sample_rank = max(rank_mu, ...
    ceil(options.random_oversampling * rank_mu));
left_rank = min(max(1, ceil(sqrt((nleft / nright) * sample_rank))), ...
    nleft);
right_rank = min(max(1, ceil(sqrt((nright / nleft) * sample_rank))), ...
    nright);

local_progress(options, 0.04, 'sampling randomized sketch');
left_sample = sample_cast(left{1}, options);
right_sample = sample_cast(right{1}, options);
left_projection = randn_like(nleft, left_rank, left_sample);
right_projection = randn_like(nright, right_rank, right_sample);
has_complex = any(cellfun(@(values) ~isreal(values), left)) || ...
    any(cellfun(@(values) ~isreal(values), right));
if has_complex
    left_projection = left_projection + 1i * ...
        randn_like(nleft, left_rank, left_sample);
    right_projection = right_projection + 1i * ...
        randn_like(nright, right_rank, right_sample);
end

products = [];
for icomponent = 1:numel(left)
    local_progress(options, 0.06 + 0.10 * (icomponent - 1) / ...
        numel(left), sprintf('sampling left projection component %d/%d', ...
        icomponent, numel(left)));
    left_values = sample_cast(left{icomponent}, options);
    compressed_left = conj(left_values) * left_projection;
    local_progress(options, 0.11 + 0.10 * (icomponent - 1) / ...
        numel(left), sprintf('sampling right projection component %d/%d', ...
        icomponent, numel(left)));
    right_values = sample_cast(right{icomponent}, options);
    compressed_right = right_values * right_projection;
    local_progress(options, 0.20 + 0.08 * (icomponent - 1) / ...
        numel(left), sprintf('sampling sketch products component %d/%d', ...
        icomponent, numel(left)));
    component_products = pair_products(compressed_left, compressed_right, ...
        @(current, total) local_product_progress( ...
        options, icomponent, numel(left), current, total));
    if isempty(products)
        products = component_products;
    else
        products = products + component_products;
    end
end
local_progress(options, 0.30, 'sampling QRCP');
ind_mu = qrcp_sample(products, rank_mu);
end

function local_product_progress(options, icomponent, ncomponents, current, total)
component_fraction = (icomponent - 1 + current / total) / ncomponents;
fraction = 0.20 + 0.08 * component_fraction;
local_progress(options, fraction, sprintf( ...
    'sampling sketch products component %d/%d: %d/%d', ...
    icomponent, ncomponents, current, total));
end

function local_progress(options, fraction, stage)
if isfield(options, 'progress') && isa(options.progress, 'function_handle')
    options.progress(fraction, stage);
end
end
