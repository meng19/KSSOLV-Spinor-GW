function coeff = product_coefficients(space, columns)
%ISDF.PRODUCT_COEFFICIENTS Return selected ISDF band-pair coefficients.
%   Coefficients are returned from the stored product_mu when available.
%   Streamed spaces reconstruct only the requested columns from the component
%   values at interpolation points.  Product columns are left-band fastest.

if nargin < 2 || isempty(columns)
    columns = 1:(space.nleft * space.nright);
end
columns = columns(:).';
npairs = space.nleft * space.nright;
if any(columns < 1 | columns > npairs | columns ~= floor(columns))
    error('ISDF:ProductColumns', 'Requested product columns are out of range.');
end
if ~isempty(space.product_mu)
    coeff = space.product_mu(:, columns);
    return;
end

nmu = numel(space.ind_mu);
block_size = min(numel(columns), max(1, floor( ...
    space.options.projection_block_elements / max(1, nmu))));
reference = space.left_mu_components{1};
coeff = complex(zeros(nmu, numel(columns), 'like', reference));
for first = 1:block_size:numel(columns)
    last = min(numel(columns), first + block_size - 1);
    block_columns = columns(first:last);
    ileft = mod(block_columns - 1, space.nleft) + 1;
    iright = floor((block_columns - 1) / space.nleft) + 1;
    block = complex(zeros(nmu, numel(block_columns), 'like', reference));
    for icomponent = 1:numel(space.left_mu_components)
        block = block + conj(space.left_mu_components{icomponent}(:, ileft)) .* ...
            space.right_mu_components{icomponent}(:, iright);
    end
    coeff(:, first:last) = block;
end
if isfield(space, 'polar_transform') && ~isempty(space.polar_transform)
    coeff = space.polar_transform * coeff;
end
end
