function [zeta_g, info] = zeta_to_g( ...
    numerator, denominator, idx_q, fftgrid, options)
%ZETA_TO_G Solve optional Gram factors and transform interpolation vectors.
%   Empty idx_q and fftgrid denote numerator values already in G space.

if isempty(idx_q) || isempty(fftgrid)
    if isempty(denominator)
        zeta_g = numerator;
        info = [];
    else
        [zeta_g, info] = stable_solve(numerator, denominator, options);
    end
    return;
end

if isempty(denominator)
    zeta_real = numerator;
    info = [];
else
    [zeta_real, info] = stable_solve(numerator, denominator, options);
end

ngrid = prod(fftgrid);
idx_q = idx_q(:);
zeta_g = complex(zeros(numel(idx_q), size(zeta_real, 2), ...
    'like', zeta_real));
for imu = 1:size(zeta_real, 2)
    zeta_grid = reshape(zeta_real(:, imu), fftgrid);
    zeta_fft = ifftn(zeta_grid);
%     zeta_fft = fftn(zeta_grid) / ngrid;
    zeta_g(:, imu) = zeta_fft(idx_q);
end
end
