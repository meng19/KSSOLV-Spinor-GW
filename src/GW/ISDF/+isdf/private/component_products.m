function [products, weight] = component_products( ...
    left, right, grid_indices, projection, sample_precision)
%COMPONENT_PRODUCTS Build or project sum_s conj(left_s).*right_s products.

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

[ngrid, nleft] = size(left{1});
[right_grid, nright] = size(right{1});
if ngrid ~= right_grid
    error('ISDF:GridMismatch', ...
        'Left and right components must share a grid.');
end
for icomponent = 1:numel(left)
    if size(left{icomponent}, 1) ~= ngrid || ...
            size(right{icomponent}, 1) ~= ngrid || ...
            size(left{icomponent}, 2) ~= nleft || ...
            size(right{icomponent}, 2) ~= nright
        error('ISDF:ComponentSizeMismatch', ...
            'All components must share grid and band dimensions.');
    end
end

if nargin < 3 || isempty(grid_indices)
    grid_indices = (1:ngrid).';
else
    grid_indices = grid_indices(:);
end
if nargin < 4
    projection = [];
end
if nargin < 5 || isempty(sample_precision)
    sample_precision = 'double';
end
if isa(left{1}, 'gpuArray') && ~isempty(projection) && ...
        ~isa(projection, 'gpuArray')
    projection = gpuArray(projection);
end

if isempty(projection)
    products = complex(zeros(numel(grid_indices), nleft * nright, ...
        'like', left{1}));
    for iright = 1:nright
        for ileft = 1:nleft
            column = ileft + (iright - 1) * nleft;
            for icomponent = 1:numel(left)
                products(:, column) = products(:, column) + ...
                    conj(left{icomponent}(grid_indices, ileft)) .* ...
                    right{icomponent}(grid_indices, iright);
            end
        end
    end
else
    if size(projection, 1) ~= nleft * nright
        error('ISDF:ProjectionSizeMismatch', ...
            ['Projection row count must match the number of ' ...
             'left-right band pairs.']);
    end
    nprojection = size(projection, 2);
    product_reference = left{1}(1);
    if strcmpi(sample_precision, 'single')
        product_reference = single(product_reference);
    end
    products = complex(zeros(numel(grid_indices), nprojection, ...
        'like', product_reference));
    ngrid_selected = numel(grid_indices);
    % Batch projection columns so that each component uses one BLAS GEMM
    % instead of one GEMM per projection.  Limit the Ngrid-by-Nright-by-Nb
    % temporary to roughly two million complex elements.
    max_batch_elements = 2e6;
    projection_block_size = min(nprojection, max(1, floor( ...
        max_batch_elements / max(1, ngrid_selected * nright))));
    for icomponent = 1:numel(left)
        left_values = left{icomponent}(grid_indices, :);
        right_values = right{icomponent}(grid_indices, :);
        if strcmpi(sample_precision, 'single')
            left_values = single(left_values);
            right_values = single(right_values);
        end
        for first_projection = 1:projection_block_size:nprojection
            last_projection = min(nprojection, ...
                first_projection + projection_block_size - 1);
            block_indices = first_projection:last_projection;
            nblock = numel(block_indices);
            projection_block = reshape(projection(:, block_indices), ...
                nleft, nright * nblock);
            transformed_left = conj(left_values) * projection_block;
            transformed_left = reshape(transformed_left, ...
                ngrid_selected, nright, nblock);
            weighted_right = transformed_left .* reshape(right_values, ...
                ngrid_selected, nright, 1);
            products(:, block_indices) = products(:, block_indices) + ...
                reshape(sum(weighted_right, 2), ngrid_selected, nblock);
        end
    end
end

if nargout > 1
    if isempty(projection)
        weight = sum(abs(products).^2, 2);
    else
        explicit = component_products(left, right, grid_indices, []);
        weight = sum(abs(explicit).^2, 2);
    end
end
end
