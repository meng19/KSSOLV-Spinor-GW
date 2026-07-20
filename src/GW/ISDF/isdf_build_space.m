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
[zeta_g, solve_info] = isdf_kernelg_current_fft(phi, psi, ind_mu, idx_q, fftgrid, options);

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
    space.product_mu = isdf_prod_states(space.phi_mu, space.psi_mu);
    return;
end

products = isdf_component_products(left_components, right_components);
space = isdf_product_space_from_products(products, idx_q, fftgrid, options);
space.left_mu_components = cell(size(left_components));
space.right_mu_components = cell(size(right_components));
for icomponent = 1:numel(left_components)
    space.left_mu_components{icomponent} = left_components{icomponent}(space.ind_mu, :);
    space.right_mu_components{icomponent} = right_components{icomponent}(space.ind_mu, :);
end
end
