function weight = component_weight(left, right)
%COMPONENT_WEIGHT Compute component-product sampling weights.

ngrid = size(left{1}, 1);
nleft = size(left{1}, 2);
nright = size(right{1}, 2);
weight = zeros(ngrid, 1, 'like', left{1});
for icomponent = 1:numel(left)
    left_i = left{icomponent};
    right_i = right{icomponent};
    if size(left_i, 1) ~= ngrid || size(right_i, 1) ~= ngrid || ...
            size(left_i, 2) ~= nleft || size(right_i, 2) ~= nright
        error('ISDF:ComponentSizeMismatch', ...
            'All components must share grid and band dimensions.');
    end
    for jcomponent = 1:numel(left)
        left_j = left{jcomponent};
        right_j = right{jcomponent};
        if size(left_j, 1) ~= ngrid || size(right_j, 1) ~= ngrid || ...
                size(left_j, 2) ~= nleft || size(right_j, 2) ~= nright
            error('ISDF:ComponentSizeMismatch', ...
                'All components must share grid and band dimensions.');
        end
        weight = weight + sum(conj(left_i) .* left_j, 2) .* ...
            sum(right_i .* conj(right_j), 2);
    end
end
weight = max(real(weight), 0);
end
