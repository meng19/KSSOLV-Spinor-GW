function space = build_space(left, right, idx_q, fftgrid, options)
%ISDF.BUILD_SPACE Build a compact component-product representation.

if nargin < 5 || isempty(options)
    options = struct();
end
if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
if numel(left) ~= numel(right)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end

ngrid = size(left{1}, 1);
nleft = size(left{1}, 2);
nright = size(right{1}, 2);
options = set_defaults(options, nleft, nright, ngrid);
if ~isfield(options, 'fftgrid') || isempty(options.fftgrid)
    options.fftgrid = fftgrid;
end

[ind_mu, products] = sample_points(left, right, options);
if isempty(products)
    product_mu = component_products(left, right, ind_mu, []);
else
    product_mu = products(ind_mu, :);
end

if strcmpi(options.sample_method, 'qrcp') && numel(left) > 1
    [zeta_real, solve_info] = stable_solve(products, product_mu, options);
    zeta_g = zeta_to_g(zeta_real, [], idx_q, fftgrid, options);
else
    [c1, c2] = product_gram(left, right, ind_mu);
    [zeta_g, solve_info] = zeta_to_g(c1, c2, idx_q, fftgrid, options);
end

space = struct();
if numel(left) == 1
    space.phi = conj(left{1});
    space.psi = right{1};
end
if strcmpi(options.sample_method, 'qrcp') && numel(left) > 1
    space.products = products;
end
space.ind_mu = ind_mu;
space.product_mu = product_mu;
space.zeta_g = zeta_g;
space.phi_mu = conj(left{1}(ind_mu, :));
space.psi_mu = right{1}(ind_mu, :);
space.left_mu_components = cell(size(left));
space.right_mu_components = cell(size(right));
for icomponent = 1:numel(left)
    space.left_mu_components{icomponent} = left{icomponent}(ind_mu, :);
    space.right_mu_components{icomponent} = right{icomponent}(ind_mu, :);
end
space.rank = numel(ind_mu);
space.options = options;
space.solve_info = solve_info;
if numel(left) > 1 && ~strcmpi(options.sample_method, 'qrcp')
    space.ngrid = ngrid;
    space.nleft = nleft;
    space.nright = nright;
end
end
