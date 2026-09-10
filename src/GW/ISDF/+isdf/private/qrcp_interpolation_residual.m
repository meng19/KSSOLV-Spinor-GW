function [residual, zeta, solve_info] = qrcp_interpolation_residual( ...
    products, ind_mu, options)
%QRCP_INTERPOLATION_RESIDUAL Relative Frobenius residual of ISDF products.

product_norm = gather_if_gpu(norm(products, 'fro'));
if product_norm == 0
    residual = 0;
    zeta = [];
    solve_info = struct();
    return;
end
product_mu = products(ind_mu, :);
[zeta, solve_info] = stable_solve(products, product_mu, options);
residual = gather_if_gpu(norm(products - zeta * product_mu, 'fro') / ...
    product_norm);
end
