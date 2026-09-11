function [products, weight] = component_products( ...
    left, right, grid_indices, projection, sample_precision, block_elements)
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
if nargin < 6 || isempty(block_elements)
    block_elements = 2e6;
end
if isa(left{1}, 'gpuArray') && ~isempty(projection) && ...
        ~isa(projection, 'gpuArray')
    projection = gpuArray(projection);
end

if isempty(projection)
    products = complex(zeros(numel(grid_indices), nleft * nright, ...
        'like', left{1}));
    % Materialize products in right-band blocks.  The former scalar loop
    % performed one MATLAB assignment per (left,right) pair.  Here each
    % block forms the same left-fastest ordering with one vectorized
    % broadcast, which is substantially faster for product_mu at selected
    % interpolation points.  Cap the temporary 3-D block near 2e6 complex
    % values; retain the scalar fallback when even one right band is larger.
    ngrid_selected = numel(grid_indices);
    right_block_size = floor(block_elements / ...
        max(1, ngrid_selected * nleft));
    if right_block_size < 1
        products = local_scalar_products(products, left, right, grid_indices);
    else
        right_block_size = min(nright, right_block_size);
        for first_right = 1:right_block_size:nright
            last_right = min(nright, first_right + right_block_size - 1);
            right_indices = first_right:last_right;
            columns = (first_right - 1) * nleft + ...
                (1:nleft * numel(right_indices));
            for icomponent = 1:numel(left)
                left_values = conj(left{icomponent}(grid_indices, :));
                right_values = right{icomponent}(grid_indices, right_indices);
                block = bsxfun(@times, reshape(left_values, ...
                    ngrid_selected, nleft, 1), reshape(right_values, ...
                    ngrid_selected, 1, numel(right_indices)));
                products(:, columns) = products(:, columns) + ...
                    reshape(block, ngrid_selected, []);
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
    % A full projected tensor has ngrid*nright*nprojection complex values.
    % Use one block whenever it fits; otherwise block only to respect the
    % temporary-memory limit.  Each block is evaluated by GEMM, with only
    % the (normally one or two) spinor components left as a loop.
    projection_block_size = min(nprojection, max(1, floor( ...
        block_elements / max(1, ngrid_selected * nright))));
    left_values = cell(1, numel(left));
    right_values = cell(1, numel(right));
    for icomponent = 1:numel(left)
        left_values{icomponent} = left{icomponent}(grid_indices, :);
        right_values{icomponent} = right{icomponent}(grid_indices, :);
        if strcmpi(sample_precision, 'single')
            left_values{icomponent} = single(left_values{icomponent});
            right_values{icomponent} = single(right_values{icomponent});
        end
    end
    for first_projection = 1:projection_block_size:nprojection
        last_projection = min(nprojection, ...
            first_projection + projection_block_size - 1);
        block_indices = first_projection:last_projection;
        nblock = numel(block_indices);
        projection_block = reshape(projection(:, block_indices), ...
            nleft, nright * nblock);
        projected_block = complex(zeros(ngrid_selected, nblock, ...
            'like', product_reference));
        for icomponent = 1:numel(left)
            transformed_left = conj(left_values{icomponent}) * projection_block;
            transformed_left = reshape(transformed_left, ...
                ngrid_selected, nright, nblock);
            projected_block = projected_block + reshape(sum( ...
                transformed_left .* reshape(right_values{icomponent}, ...
                ngrid_selected, nright, 1), 2), ngrid_selected, nblock);
        end
        products(:, block_indices) = projected_block;
    end
end

if nargout > 1
    if isempty(projection)
        weight = sum(abs(products).^2, 2);
    else
        explicit = feval(mfilename, left, right, grid_indices, [], ...
            sample_precision, block_elements);
        weight = sum(abs(explicit).^2, 2);
    end
end
end

function products = local_scalar_products(products, left, right, grid_indices)
% Keep peak temporary memory bounded for a full-grid explicit QRCP matrix.
for iright = 1:size(right{1}, 2)
    for ileft = 1:size(left{1}, 2)
        column = ileft + (iright - 1) * size(left{1}, 2);
        for icomponent = 1:numel(left)
            products(:, column) = products(:, column) + ...
                conj(left{icomponent}(grid_indices, ileft)) .* ...
                right{icomponent}(grid_indices, iright);
        end
    end
end
end
