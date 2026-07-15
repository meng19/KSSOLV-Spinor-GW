function weight = isdf_product_weight(phi, psi, options)
%ISDF_PRODUCT_WEIGHT Build grid weights for interpolation point sampling.

rho_phi = sum(abs(phi).^2, 2);
rho_psi = sum(abs(psi).^2, 2);

switch lower(options.weight)
    case {'prod', 'multiply'}
        weight = rho_phi .* rho_psi;
    case 'add'
        weight = rho_phi + rho_psi;
    case 'power'
        weight = (rho_phi .* rho_psi).^(options.power / 2);
    otherwise
        error('ISDF:UnknownWeight', ...
            'Unknown ISDF weight "%s". Supported weights: prod, add, power.', options.weight);
end

weight = real(weight(:));
weight(~isfinite(weight) | weight < 0) = 0;
if all(weight == 0)
    weight = ones(size(weight));
end
end
