function products = isdf_component_products(left_components, right_components)
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

products = zeros(ngrid, nleft * nright);
for iright = 1:nright
    for ileft = 1:nleft
        col = ileft + (iright - 1) * nleft;
        product = zeros(ngrid, 1);
        for icomponent = 1:ncomponents
            if size(left_components{icomponent}, 1) ~= ngrid || ...
                    size(right_components{icomponent}, 1) ~= ngrid
                error('ISDF:GridMismatch', ...
                    'All components must share a grid.');
            end
            product = product + conj(left_components{icomponent}(:, ileft)) .* ...
                right_components{icomponent}(:, iright);
        end
        products(:, col) = product;
    end
end
end
