function space = isdf_build_space(phi, psi, idx_q, fftgrid, options)
%ISDF_BUILD_SPACE Build a compact ISDF product-space representation.

if nargin < 5 || isempty(options)
    options = struct();
end

if iscell(phi) || iscell(psi)
    space = isdf_build_component_space(phi, psi, idx_q, fftgrid, options);
    return;
end

options = isdf_set_defaults(options, size(phi, 2), size(psi, 2), size(phi, 1));
if ~isfield(options, 'fftgrid') || isempty(options.fftgrid)
    options.fftgrid = fftgrid;
end
ind_mu = isdf_indices(phi, psi, options);
[c1, c2] = isdf_product_gram({conj(phi)}, {psi}, ind_mu);
[zeta_g, solve_info] = isdf_zeta_g_from_product_gram(c1, c2, idx_q, fftgrid, options);

space = struct();
space.phi = phi;
space.psi = psi;
space.ind_mu = ind_mu;
space.zeta_g = zeta_g;
space.phi_mu = phi(ind_mu, :);
space.psi_mu = psi(ind_mu, :);
space.rank = numel(ind_mu);
space.options = options;
space.solve_info = solve_info;
end

function space = isdf_build_component_space(left_components, right_components, idx_q, fftgrid, options)
if ~iscell(left_components)
    left_components = {left_components};
end
if ~iscell(right_components)
    right_components = {right_components};
end
if numel(left_components) ~= numel(right_components)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end

if numel(left_components) == 1
    phi = conj(left_components{1});
    psi = right_components{1};
    space = isdf_build_space(phi, psi, idx_q, fftgrid, options);
    space.left_mu_components = {left_components{1}(space.ind_mu, :)};
    space.right_mu_components = {space.psi_mu};
    space.product_mu = isdf_component_products( ...
        space.left_mu_components, space.right_mu_components);
    return;
end

ngrid = size(left_components{1}, 1);
nleft = size(left_components{1}, 2);
nright = size(right_components{1}, 2);
options = isdf_set_defaults(options, nleft, nright, ngrid);
if ~isfield(options, 'fftgrid') || isempty(options.fftgrid)
    options.fftgrid = fftgrid;
end

switch lower(options.sample_method)
    case 'qrcp'
        products = isdf_component_products(left_components, right_components);
        space = isdf_product_space_from_products(products, idx_q, fftgrid, options);
    otherwise
        space = isdf_build_component_space_matrix_free( ...
            left_components, right_components, idx_q, fftgrid, options);
end
space.left_mu_components = cell(size(left_components));
space.right_mu_components = cell(size(right_components));
for icomponent = 1:numel(left_components)
    space.left_mu_components{icomponent} = left_components{icomponent}(space.ind_mu, :);
    space.right_mu_components{icomponent} = right_components{icomponent}(space.ind_mu, :);
end
end

function space = isdf_build_component_space_matrix_free(left_components, right_components, idx_q, fftgrid, options)
ngrid = size(left_components{1}, 1);
nleft = size(left_components{1}, 2);
nright = size(right_components{1}, 2);

rng(options.seed, 'twister');
switch lower(options.sample_method)
    case {'qrcp_randomized', 'randomized_qrcp', 'default'}
        npairs = nleft * nright;
        projection_rank = max(options.rank, ceil(options.random_oversampling * options.rank));
        projection_rank = min(projection_rank, npairs);
        random_projection = randn(npairs, projection_rank);
        has_complex = any(cellfun(@(component) ~isreal(component), left_components)) || ...
            any(cellfun(@(component) ~isreal(component), right_components));
        if has_complex
            random_projection = random_projection + 1i * randn(npairs, projection_rank);
        end
        compressed_products = isdf_component_products_projected( ...
            left_components, right_components, random_projection);
        ind_mu = isdf_qrcp_indices_from_products(compressed_products, options.rank);
    case 'kmeans'
        weight = isdf_component_product_weight(left_components, right_components);
        ind_mu = isdf_kmeans_indices(weight, options);
    otherwise
        error('ISDF:UnknownSampleMethod', ...
            ['Unknown ISDF sample_method "%s". Supported methods: qrcp, ' ...
            'qrcp_randomized, kmeans, default.'], options.sample_method);
end

product_mu = isdf_component_products(left_components, right_components, ind_mu);
[c1, c2] = isdf_product_gram(left_components, right_components, ind_mu);
[zeta_g, solve_info] = isdf_zeta_g_from_product_gram(c1, c2, idx_q, fftgrid, options);

space = struct();
space.ind_mu = ind_mu;
space.product_mu = product_mu;
space.zeta_g = zeta_g;
space.rank = numel(ind_mu);
space.options = options;
space.solve_info = solve_info;
space.ngrid = ngrid;
space.nleft = nleft;
space.nright = nright;
end
