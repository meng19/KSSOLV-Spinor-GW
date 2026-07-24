function space = isdf_product_space_from_products(products, idx_q, fftgrid, options)
%ISDF_PRODUCT_SPACE_FROM_PRODUCTS Build ISDF space from preformed products.

if nargin < 4 || isempty(options)
    options = struct();
end

[ngrid, npairs] = size(products);
options = isdf_set_defaults(options, 1, npairs, ngrid);
if ~isfield(options, 'fftgrid') || isempty(options.fftgrid)
    options.fftgrid = fftgrid;
end

rng(options.seed, 'twister');
switch lower(options.sample_method)
    case 'qrcp'
        ind_mu = isdf_qrcp_indices_from_products(products, options.rank);
    case {'qrcp_randomized', 'randomized_qrcp', 'default'}
        projection_rank = max(options.rank, ceil(options.random_oversampling * options.rank));
        projection_rank = min(projection_rank, npairs);
        random_projection = randn(npairs, projection_rank);
        if ~isreal(products)
            random_projection = random_projection + 1i * randn(npairs, projection_rank);
        end
        compressed_products = products * random_projection;
        ind_mu = isdf_qrcp_indices_from_products(compressed_products, options.rank);
    case 'kmeans'
        weight = sum(abs(products).^2, 2);
        ind_mu = isdf_kmeans_indices(weight, options);
    otherwise
        error('ISDF:UnknownSampleMethod', ...
            ['Unknown ISDF sample_method "%s". Supported methods: qrcp, ' ...
            'qrcp_randomized, kmeans, default.'], options.sample_method);
end

product_mu = products(ind_mu, :);
[zeta_real, solve_info] = isdf_stable_right_solve(products, product_mu, options);
zeta_g = isdf_zeta_real_to_g(zeta_real, idx_q, fftgrid);

space = struct();
space.products = products;
space.ind_mu = ind_mu;
space.product_mu = product_mu;
space.zeta_g = zeta_g;
space.rank = numel(ind_mu);
space.options = options;
space.solve_info = solve_info;
end
