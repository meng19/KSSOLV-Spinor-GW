function ind_mu = isdf_indices(phi, psi, options)
%ISDF_INDICES Select interpolation points for product states.

options = isdf_set_defaults(options, size(phi, 2), size(psi, 2), size(phi, 1));
rng(options.seed, 'twister');

switch lower(options.sample_method)
    case 'qrcp'
        products = isdf_prod_states(phi, psi);
        [~, ~, pivots] = qr(products', 0);
        ind_mu = pivots(1:options.rank);
        ind_mu = reshape(ind_mu, 1, []);
    otherwise
        error('ISDF:UnknownSampleMethod', ...
            'Unknown ISDF sample_method "%s". Supported method: qrcp.', options.sample_method);
end
end
