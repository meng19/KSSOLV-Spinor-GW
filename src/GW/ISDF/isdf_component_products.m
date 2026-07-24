function products = isdf_component_products(left_components, right_components, grid_indices)
%ISDF_COMPONENT_PRODUCTS Build summed component pair products.

if ~iscell(left_components)
    left_components = {left_components};
end
if ~iscell(right_components)
    right_components = {right_components};
end

ncomponents = numel(left_components);
if ncomponents ~= numel(right_components)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end

[ngrid, nleft] = size(left_components{1});
[ngrid_right, nright] = size(right_components{1});
if ngrid ~= ngrid_right
    error('ISDF:GridMismatch', 'Left and right components must share a grid.');
end
if nargin < 3 || isempty(grid_indices)
    grid_indices = (1:ngrid).';
else
    grid_indices = grid_indices(:);
end

products = zeros(numel(grid_indices), nleft * nright);
for iright = 1:nright
    for ileft = 1:nleft
        col = ileft + (iright - 1) * nleft;
        product = zeros(numel(grid_indices), 1);
        for icomponent = 1:ncomponents
            if size(left_components{icomponent}, 1) ~= ngrid || ...
                    size(right_components{icomponent}, 1) ~= ngrid || ...
                    size(left_components{icomponent}, 2) ~= nleft || ...
                    size(right_components{icomponent}, 2) ~= nright
                error('ISDF:ComponentSizeMismatch', ...
                    'All components must share grid and band dimensions.');
            end
            left_values = left_components{icomponent}(grid_indices, ileft);
            right_values = right_components{icomponent}(grid_indices, iright);
            product = product + conj(left_values) .* right_values;
        end
        products(:, col) = product;
    end
end
end
