function gme = isdf_epsilon_matrix_elements_from_real(phi, psi, idx_q, fftgrid, options)
%ISDF_EPSILON_MATRIX_ELEMENTS_FROM_REAL Approximate epsilon matrix elements.
%   phi and psi are real-space wavefunction factors on the FFT grid. For
%   epsilon, pass phi = conj(valence_real) and psi = conduction_real.

options = isdf_set_defaults(options, size(phi, 2), size(psi, 2), size(phi, 1));

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
zetag = isdf_kernelg_current_fft(phi, psi, ind_mu, idx_q, fftgrid);

rank_mu = length(ind_mu);
gme = zeros(length(idx_q), nphi, npsi);
for ipsi = 1:npsi
    for iphi = 1:nphi
        coeff = phi(ind_mu, iphi) .* psi(ind_mu, ipsi);
        gme(:, iphi, ipsi) = zetag * reshape(coeff, rank_mu, 1);
    end
end
end
