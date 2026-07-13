function products = isdf_prod_states(phi, psi)
%ISDF_PROD_STATES Build pair-product states phi_i(r) * psi_j(r).

[nphi_grid, nphi] = size(phi);
[npsi_grid, npsi] = size(psi);
if nphi_grid ~= npsi_grid
    error('ISDF:GridMismatch', 'phi and psi must have the same number of grid points.');
end

products = zeros(nphi_grid, nphi * npsi);
for ipsi = 1:npsi
    for iphi = 1:nphi
        products(:, iphi + (ipsi - 1) * nphi) = phi(:, iphi) .* psi(:, ipsi);
    end
end
end
