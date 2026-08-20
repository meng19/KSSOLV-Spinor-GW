function weight = scalar_weight(left, right, options)
%SCALAR_WEIGHT Compute scalar wavefunction sampling weights.

rho_left = sum(abs(left).^2, 2);
rho_right = sum(abs(right).^2, 2);
switch lower(options.weight)
    case {'prod', 'multiply'}
        weight = rho_left .* rho_right;
    case 'add'
        weight = rho_left + rho_right;
    case 'power'
        weight = (rho_left .* rho_right).^(options.power / 2);
    otherwise
        error('ISDF:UnknownWeight', ...
            ['Unknown ISDF weight "%s". Supported weights: ' ...
             'prod, add, power.'], options.weight);
end
weight = real(weight(:));
weight = gather_if_gpu(weight);
weight(~isfinite(weight) | weight < 0) = 0;
if all(weight == 0)
    weight = ones(size(weight));
end
end
