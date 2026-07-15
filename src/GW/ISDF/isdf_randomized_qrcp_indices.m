function ind_mu = isdf_randomized_qrcp_indices(phi, psi, options)
%ISDF_RANDOMIZED_QRCP_INDICES Select pivots after random product compression.

[~, nphi] = size(phi);
[~, npsi] = size(psi);

rank_mu = options.rank;
rsamp = max(rank_mu, ceil(options.random_oversampling * rank_mu));
rphi = min(max(1, ceil(sqrt((nphi / npsi) * rsamp))), nphi);
rpsi = min(max(1, ceil(sqrt((npsi / nphi) * rsamp))), npsi);

gphi = randn(nphi, rphi);
gpsi = randn(npsi, rpsi);
if ~isreal(phi) || ~isreal(psi)
    gphi = gphi + 1i * randn(nphi, rphi);
    gpsi = gpsi + 1i * randn(npsi, rpsi);
end

phi_g = phi * gphi;
psi_g = psi * gpsi;
compressed_products = isdf_prod_states(phi_g, psi_g);
ind_mu = isdf_qrcp_indices_from_products(compressed_products, rank_mu);
end
