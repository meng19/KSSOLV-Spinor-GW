function zeta_g = isdf_zeta_real_to_g(zeta_real, idx_q, fftgrid)
%ISDF_ZETA_REAL_TO_G Transform real-space ISDF interpolation vectors to G space.

ngrid = prod(fftgrid);
idx_q = idx_q(:);
zeta_g = zeros(length(idx_q), size(zeta_real, 2));
for imu = 1:size(zeta_real, 2)
    zeta_grid = reshape(zeta_real(:, imu), fftgrid);
    zeta_fft = fftn(zeta_grid) / ngrid;
    zeta_g(:, imu) = zeta_fft(idx_q);
end
end
