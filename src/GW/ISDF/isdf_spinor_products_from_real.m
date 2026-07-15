function products = isdf_spinor_products_from_real(left_components, right_components)
%ISDF_SPINOR_PRODUCTS_FROM_REAL Build summed spinor pair products.

nspinor = numel(left_components);
if nspinor ~= numel(right_components)
    error('ISDF:SpinorComponentMismatch', ...
        'Left and right spinor component counts must match.');
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
        for ispinor = 1:nspinor
            if size(left_components{ispinor}, 1) ~= ngrid || ...
                    size(right_components{ispinor}, 1) ~= ngrid
                error('ISDF:GridMismatch', ...
                    'All spinor components must share a grid.');
            end
            product = product + conj(left_components{ispinor}(:, ileft)) .* ...
                right_components{ispinor}(:, iright);
        end
        products(:, col) = product;
    end
end
end
