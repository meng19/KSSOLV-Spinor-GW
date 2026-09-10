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
options.weight_was_set = isfield(options, 'weight') && ...
    ~isempty(options.weight);
options = set_defaults(options, nleft, nright, ngrid);
if ~isfield(options, 'fftgrid') || isempty(options.fftgrid)
    options.fftgrid = fftgrid;
end
local_progress(options, 0.02, 'setup');
[ind_mu, products, adaptive_state] = sample_points(left, right, options);
local_progress(options, 0.35, 'sample points');
adaptive_info = struct('enabled', false, 'residual', NaN, ...
    'initial_rank', options.rank, 'final_rank', options.rank, ...
    'tolerance', NaN, 'reached_tolerance', false);
if adaptive_state.enabled
    adaptive_info = rmfield(adaptive_state, {'has_zeta', 'zeta', 'solve_info'});
    options.rank = adaptive_info.final_rank;
    options.rank_source = 'adaptive_qrcp';
end
local_print_rank(options, ngrid, nleft, nright);
if isempty(products)
    product_mu = component_products(left, right, ind_mu, []);
else
    product_mu = products(ind_mu, :);
end
local_progress(options, 0.45, 'product values');

if adaptive_state.enabled && adaptive_state.has_zeta
    zeta_real = adaptive_state.zeta;
    solve_info = adaptive_state.solve_info;
    zeta_g = zeta_to_g(zeta_real, [], idx_q, fftgrid, options);
elseif strcmpi(options.sample_method, 'qrcp') && numel(left) > 1
    local_progress(options, 0.50, 'interpolation solve');
    [zeta_real, solve_info] = stable_solve(products, product_mu, options);
    local_progress(options, 0.82, 'Fourier transform');
    zeta_g = zeta_to_g(zeta_real, [], idx_q, fftgrid, options);
else
    local_progress(options, 0.50, 'product Gram');
    [c1, c2] = product_gram(left, right, ind_mu);
    local_progress(options, 0.78, 'interpolation solve and FFT');
    [zeta_g, solve_info] = zeta_to_g(c1, c2, idx_q, fftgrid, options);
end
local_progress(options, 0.95, 'complete');

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
space.adaptive_info = adaptive_info;
if numel(left) > 1 && ~strcmpi(options.sample_method, 'qrcp')
    space.ngrid = ngrid;
    space.nleft = nleft;
    space.nright = nright;
end
end

function local_progress(options, fraction, stage)
if isfield(options, 'progress') && isa(options.progress, 'function_handle')
    options.progress(fraction, stage);
end
end

function local_print_rank(options, ngrid, nleft, nright)
if isfield(options, 'print_rank') && ~options.print_rank
    return;
end

persistent printed_keys;
if isempty(printed_keys)
    printed_keys = {};
end

key = sprintf('%s:%s:%d:%d:%d:%d:%d:%d:%.16g', ...
    lower(char(options.sample_method)), lower(char(options.rank_source)), ...
    ngrid, nleft, nright, options.rank, options.recommended_rank, ...
    options.max_rank, options.rank_ratio);
if any(strcmp(printed_keys, key))
    return;
end
printed_keys{end + 1} = key;

fprintf(['\nISDF rank: rank = %d, recommended = ' ...
    'ceil(sqrt(%d*%d)*%.3g) = %d\n'], ...
    options.rank, nleft, nright, options.rank_ratio, ...
    options.recommended_rank);
end
