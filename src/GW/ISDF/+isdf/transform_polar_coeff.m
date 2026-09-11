function coeff = transform_polar_coeff(coeff, transform)
%TRANSFORM_POLAR_COEFF Change reduced polarizability to a whitened basis.

if isempty(transform)
    return;
end
for ifreq = 1:size(coeff, 3)
    coeff(:, :, ifreq) = transform * coeff(:, :, ifreq) * transform';
end
end
