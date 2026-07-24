function gme = isdf_matrix_elements_from_real(phi, psi, idx_q, fftgrid, options)
%ISDF_MATRIX_ELEMENTS_FROM_REAL Approximate FFTs of real-space products.
%   phi and psi are real-space wavefunction factors on the FFT grid. The
%   returned array stores FFT_G(phi(:, i) .* psi(:, j)) for all pairs.

options = isdf_set_defaults(options, size(phi, 2), size(psi, 2), size(phi, 1));
if ~isfield(options, 'fftgrid') || isempty(options.fftgrid)
    options.fftgrid = fftgrid;
end

[nphi_grid, nphi] = size(phi);
[npsi_grid, npsi] = size(psi);
if nphi_grid ~= npsi_grid
    error('ISDF:GridMismatch', 'phi and psi must have the same number of grid points.');
end
if prod(fftgrid) ~= nphi_grid
    error('ISDF:GridSizeMismatch', 'prod(fftgrid) must match the number of rows in phi and psi.');
end

idx_q = idx_q(:);
ind_mu = isdf_indices(phi, psi, options);
[c1, c2] = isdf_product_gram({conj(phi)}, {psi}, ind_mu);
zetag = isdf_zeta_g_from_product_gram(c1, c2, idx_q, fftgrid, options);

rank_mu = length(ind_mu);
gme = zeros(length(idx_q), nphi, npsi);
for ipsi = 1:npsi
    for iphi = 1:nphi
        coeff = phi(ind_mu, iphi) .* psi(ind_mu, ipsi);
        gme(:, iphi, ipsi) = zetag * reshape(coeff, rank_mu, 1);
    end
end
end
