function weight = isdf_component_product_weight(left_components, right_components)
%ISDF_COMPONENT_PRODUCT_WEIGHT Compute row norms of component product states.
%   This returns sum(abs(P).^2, 2) for
%
%       P(r,ij) = sum_s conj(left_s(r,i)) * right_s(r,j)
%
%   without explicitly forming the full product matrix P.

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

ngrid = size(left_components{1}, 1);
weight = zeros(ngrid, 1);
for icomponent = 1:numel(left_components)
    left_i = left_components{icomponent};
    right_i = right_components{icomponent};
    if size(left_i, 1) ~= ngrid || size(right_i, 1) ~= ngrid || ...
            size(left_i, 2) ~= size(left_components{1}, 2) || ...
            size(right_i, 2) ~= size(right_components{1}, 2)
        error('ISDF:ComponentSizeMismatch', ...
            'All components must share grid and band dimensions.');
    end
    for jcomponent = 1:numel(left_components)
        left_j = left_components{jcomponent};
        right_j = right_components{jcomponent};
        if size(left_j, 1) ~= ngrid || size(right_j, 1) ~= ngrid || ...
                size(left_j, 2) ~= size(left_components{1}, 2) || ...
                size(right_j, 2) ~= size(right_components{1}, 2)
            error('ISDF:ComponentSizeMismatch', ...
                'All components must share grid and band dimensions.');
        end
        weight = weight + ...
            sum(conj(left_i) .* left_j, 2) .* ...
            sum(right_i .* conj(right_j), 2);
    end
end
weight = max(real(weight), 0);
end
