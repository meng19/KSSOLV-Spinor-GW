function [products, weight] = component_products( ...
    left, right, grid_indices, projection)
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

if isempty(projection)
    products = zeros(numel(grid_indices), nleft * nright);
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
    products = zeros(numel(grid_indices), size(projection, 2));
    for iprojection = 1:size(projection, 2)
        projection_matrix = reshape(projection(:, iprojection), ...
            nleft, nright);
        for icomponent = 1:numel(left)
            left_values = left{icomponent}(grid_indices, :);
            right_values = right{icomponent}(grid_indices, :);
            products(:, iprojection) = products(:, iprojection) + ...
                sum((conj(left_values) * projection_matrix) .* ...
                right_values, 2);
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
