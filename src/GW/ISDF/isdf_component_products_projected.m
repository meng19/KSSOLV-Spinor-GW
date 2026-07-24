function projected_products = isdf_component_products_projected(left_components, right_components, projection)
%ISDF_COMPONENT_PRODUCTS_PROJECTED Compute P * projection without forming P.
%   P is the component product matrix
%
%       P(r,ij) = sum_s conj(left_s(r,i)) * right_s(r,j).

if ~iscell(left_components)
    left_components = {left_components};
end
if ~iscell(right_components)
    right_components = {right_components};
end
if numel(left_components) ~= numel(right_components)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end

[ngrid, nleft] = size(left_components{1});
[ngrid_right, nright] = size(right_components{1});
if ngrid ~= ngrid_right
    error('ISDF:GridMismatch', 'Left and right components must share a grid.');
end
if size(projection, 1) ~= nleft * nright
    error('ISDF:ProjectionSizeMismatch', ...
        'Projection row count must match the number of left-right band pairs.');
end

projected_products = zeros(ngrid, size(projection, 2));
for iproj = 1:size(projection, 2)
    projection_matrix = reshape(projection(:, iproj), nleft, nright);
    projected_product = zeros(ngrid, 1);
    for icomponent = 1:numel(left_components)
        left = left_components{icomponent};
        right = right_components{icomponent};
        if size(left, 1) ~= ngrid || size(right, 1) ~= ngrid || ...
                size(left, 2) ~= nleft || size(right, 2) ~= nright
            error('ISDF:ComponentSizeMismatch', ...
                'All components must share grid and band dimensions.');
        end
        projected_product = projected_product + ...
            sum((conj(left) * projection_matrix) .* right, 2);
    end
    projected_products(:, iproj) = projected_product;
end
end
