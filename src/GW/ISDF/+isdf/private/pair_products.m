function products = pair_products(left, right)
%PAIR_PRODUCTS Materialize scalar pair products by grid row.

[ngrid, nleft] = size(left);
nright = size(right, 2);
products = complex(zeros(ngrid, nleft * nright, 'like', left));
for iright = 1:nright
    for ileft = 1:nleft
        products(:, ileft + (iright - 1) * nleft) = ...
            left(:, ileft) .* right(:, iright);
    end
end
end
