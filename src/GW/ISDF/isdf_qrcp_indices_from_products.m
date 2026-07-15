function ind_mu = isdf_qrcp_indices_from_products(products, rank_mu)
%ISDF_QRCP_INDICES_FROM_PRODUCTS Select grid pivots from product states.

[~, ~, pivots] = qr(products', 0);
ind_mu = reshape(pivots(1:rank_mu), 1, []);
end
