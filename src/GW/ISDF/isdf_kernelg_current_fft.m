function zetag = isdf_kernelg_current_fft(phi, psi, ind_mu, idx_q, fftgrid)
%ISDF_KERNELG_CURRENT_FFT Build ISDF helper functions in current FFT convention.

[ngrid, ~] = size(phi);
rank_mu = length(ind_mu);
idx_q = idx_q(:);

phi_mu = phi(ind_mu, :);
psi_mu = psi(ind_mu, :);
c2 = (phi_mu * phi_mu') .* (psi_mu * psi_mu');
c1 = (phi * phi_mu') .* (psi * psi_mu');
zeta_real = c1 / c2;

zetag = zeros(length(idx_q), rank_mu);
for imu = 1:rank_mu
    zeta_grid = reshape(zeta_real(:, imu), fftgrid);
    zeta_fft = fftn(zeta_grid) / ngrid;
    zetag(:, imu) = zeta_fft(idx_q);
end
end
