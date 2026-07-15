function ind_mu = isdf_indices(phi, psi, options)
%ISDF_INDICES Select interpolation points for product states.

options = isdf_set_defaults(options, size(phi, 2), size(psi, 2), size(phi, 1));
rng(options.seed, 'twister');

switch lower(options.sample_method)
    case 'qrcp'
        products = isdf_prod_states(phi, psi);
        ind_mu = isdf_qrcp_indices_from_products(products, options.rank);
    case {'qrcp_randomized', 'randomized_qrcp'}
        ind_mu = isdf_randomized_qrcp_indices(phi, psi, options);
    case 'kmeans'
        weight = isdf_product_weight(phi, psi, options);
        ind_mu = isdf_kmeans_indices(weight, options);
    case 'default'
        ind_mu = isdf_randomized_qrcp_indices(phi, psi, options);
    otherwise
        error('ISDF:UnknownSampleMethod', ...
            ['Unknown ISDF sample_method "%s". Supported methods: qrcp, ' ...
            'qrcp_randomized, kmeans, default.'], options.sample_method);
end
end
