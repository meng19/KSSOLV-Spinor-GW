function space = build_space(left, right, idx_q, fftgrid, options)
%ISDF.BUILD_SPACE Build a compact component-product representation.

if nargin < 5 || isempty(options)
    options = struct();
end
if ~iscell(left)
    left = {left};
end
if isfield(options, 'right_component_fun') && ...
        ~isempty(options.right_component_fun)
    space = local_build_space_streamed(left, idx_q, fftgrid, options);
    return;
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
options = set_defaults(options, nleft, nright, ngrid, fftgrid);
local_print_rank(options, ngrid, nleft, nright);
local_progress(options, 0.02, 'setup');
if options.swap_left_right && options.adaptive_rank_enable
    error('ISDF:SwapAdaptiveRank', ...
        'swap_left_right is not supported with adaptive_rank_enable.');
end
if options.swap_left_right
    % Keep the physical representation in the original left-right order.
    % Only the QRCP/k-means point-selection control sees the conjugate
    % right-left product space.
    [ind_mu, ~, adaptive_state] = sample_points(right, left, options);
    products = [];
else
    [ind_mu, products, adaptive_state] = sample_points(left, right, options);
end
local_progress(options, 0.35, 'sample points');
adaptive_info = struct('enabled', false, 'residual', NaN, ...
    'initial_rank', options.rank, 'final_rank', options.rank, ...
    'tolerance', NaN, 'reached_tolerance', false);
if adaptive_state.enabled
    adaptive_info = rmfield(adaptive_state, {'has_zeta', 'zeta', 'solve_info'});
    options.rank = adaptive_info.final_rank;
    options.rank_source = 'adaptive_qrcp';
end
materialize_product_mu = local_materialize_product_mu(options, nleft, nright);
if isempty(products) && materialize_product_mu
    product_mu = component_products(left, right, ind_mu, [], ...
        options.sample_precision, options.projection_block_elements);
elseif ~isempty(products)
    product_mu = products(ind_mu, :);
else
    product_mu = [];
end
local_progress(options, 0.45, 'product values');

transform = [];
if adaptive_state.enabled && adaptive_state.has_zeta
    if strcmp(options.interpolation_solver, 'svd_whiten')
        error('ISDF:SVDWhitenAdaptive', ...
            'svd_whiten is not supported with adaptive_rank_enable.');
    end
    zeta_real = adaptive_state.zeta;
    solve_info = adaptive_state.solve_info;
    zeta_g = zeta_to_g(zeta_real, [], idx_q, fftgrid, options);
elseif strcmpi(options.sample_method, 'qrcp') && numel(left) > 1 && ...
        ~isempty(products)
    if strcmp(options.interpolation_solver, 'svd_whiten')
        local_progress(options, 0.50, 'product Gram and SVD whitening');
        [c1, c2] = product_gram(left, right, ind_mu);
        local_progress(options, 0.78, 'whitened interpolation and FFT');
        [zeta_g, product_mu, transform, solve_info] = svd_whiten( ...
            c1, c2, product_mu, idx_q, fftgrid, options);
    else
        local_progress(options, 0.50, 'interpolation solve');
        [zeta_real, solve_info] = stable_solve(products, product_mu, options);
        local_progress(options, 0.82, 'Fourier transform');
        zeta_g = zeta_to_g(zeta_real, [], idx_q, fftgrid, options);
    end
else
    local_progress(options, 0.50, 'product Gram');
    if strcmp(options.interpolation_solver, 'svd_whiten')
        [c1, c2] = product_gram(left, right, ind_mu);
        local_progress(options, 0.78, 'SVD whitening and FFT');
        [zeta_g, product_mu, transform, solve_info] = svd_whiten( ...
            c1, c2, product_mu, idx_q, fftgrid, options);
%         [zeta_g, ~, transform, solve_info] = svd_whiten( ...
%             c1, c2, product_mu, idx_q, fftgrid, options);
    else
        local_progress(options, 0.78, 'interpolation solve and FFT');
        if local_stream_product_gram(options, ngrid)
            [c1_g, c2] = product_gram(left, right, ind_mu, ...
                idx_q, fftgrid, options);
            [zeta_g, solve_info] = zeta_to_g(c1_g, c2, [], [], options);
        else
            [c1, c2] = product_gram(left, right, ind_mu);
            [zeta_g, solve_info] = zeta_to_g(c1, c2, ...
                idx_q, fftgrid, options);
        end
    end
end

% Exact QRCP necessarily materializes products while selecting pivots, but
% the final space can still release its selected coefficient matrix when
% streamed storage was requested.  product_coefficients then reconstructs it.
if ~materialize_product_mu
    product_mu = [];
end

space = struct();
if numel(left) == 1
    space.phi = conj(left{1});
    space.psi = right{1};
end
if strcmpi(options.sample_method, 'qrcp') && numel(left) > 1 && ...
        ~isempty(products)
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
space.sampling_rank = numel(ind_mu);
if ~isempty(transform)
    space.rank = size(transform, 1);
    space.polar_transform = transform;
    fprintf('ISDF SVD whitening: retained %d/%d modes (cutoff %.3g, ratio %.3g)\n', ...
        space.rank, space.sampling_rank, options.svd_cutoff, options.svd_ratio);
end
space.options = options;
space.solve_info = solve_info;
space.adaptive_info = adaptive_info;
space.sampling_swapped_left_right = options.swap_left_right;
space.ngrid = ngrid;
space.nleft = nleft;
space.nright = nright;
local_progress(options, 1, 'complete');
end

function space = local_build_space_streamed(left, idx_q, fftgrid, options)
% Build one global space while generating right-band real components in
% bounded blocks.  Only sampled right components are retained afterwards.
provider = options.right_component_fun;
if ~isfield(options, 'right_bands') || isempty(options.right_bands)
    error('ISDF:StreamedRightBands', ...
        'right_component_fun requires right_bands.');
end
right_bands = options.right_bands(:).';
ngrid = size(left{1}, 1);
nleft = size(left{1}, 2);
nright = numel(right_bands);
options.weight_was_set = isfield(options, 'weight') && ...
    ~isempty(options.weight);
options = set_defaults(options, nleft, nright, ngrid, fftgrid);
local_print_rank(options, ngrid, nleft, nright);
if ~strcmpi(options.sample_method, 'kmeans')
    error('ISDF:StreamedRightSampling', ...
        'Streamed right components currently require kmeans sampling.');
end
if options.swap_left_right || options.adaptive_rank_enable
    error('ISDF:StreamedRightControls', ...
        ['Streamed right components do not support swap_left_right ' ...
        'or adaptive_rank_enable.']);
end
if ~strcmp(options.interpolation_solver, 'direct')
    error('ISDF:StreamedRightSolver', ...
        'Streamed right components require interpolation_solver=direct.');
end
if ~strcmp(options.product_mu_storage, 'stream') || ...
        ~local_stream_product_gram(options, ngrid)
    error('ISDF:StreamedRightStorage', ...
        ['Streamed right components require product_mu_storage=stream ' ...
        'and streamed product Gram storage.']);
end
if ~isfinite(options.right_band_block_size)
    error('ISDF:StreamedRightBlockSize', ...
        'right_band_block_size must be finite for streamed right components.');
end
block_size = min(nright, max(1, ceil(options.right_band_block_size)));
local_progress(options, 0.02, 'setup');
rng(options.seed, 'twister');
weight = local_streamed_weight( ...
    left, provider, right_bands, block_size, options);
ind_mu = kmeans_sample(weight, options);
local_progress(options, 0.35, 'sample points');

nmu = numel(ind_mu);
right_mu = cell(1, numel(left));
for icomponent = 1:numel(left)
    right_mu{icomponent} = complex(zeros(nmu, nright, 'like', left{1}));
end
c1_g = [];
c2 = [];
for first = 1:block_size:nright
    last = min(nright, first + block_size - 1);
    bands = right_bands(first:last);
    right = provider(bands);
    if ~iscell(right)
        right = {right};
    end
    [c1_block, c2_block] = product_gram( ...
        left, right, ind_mu, idx_q, fftgrid, options);
    if isempty(c1_g)
        c1_g = c1_block;
        c2 = c2_block;
    else
        c1_g = c1_g + c1_block;
        c2 = c2 + c2_block;
    end
    for icomponent = 1:numel(left)
        right_mu{icomponent}(:, first:last) = right{icomponent}(ind_mu, :);
    end
    local_progress(options, 0.35 + 0.43 * last / nright, ...
        sprintf('streamed product Gram bands %d/%d', last, nright));
end
local_progress(options, 0.78, 'interpolation solve and FFT');
[zeta_g, solve_info] = zeta_to_g(c1_g, c2, [], [], options);

space = struct();
space.ind_mu = ind_mu;
space.product_mu = [];
space.zeta_g = zeta_g;
space.phi_mu = conj(left{1}(ind_mu, :));
space.psi_mu = right_mu{1};
space.left_mu_components = cell(size(left));
space.right_mu_components = right_mu;
for icomponent = 1:numel(left)
    space.left_mu_components{icomponent} = left{icomponent}(ind_mu, :);
end
space.rank = nmu;
space.sampling_rank = nmu;
if isfield(options, 'right_component_fun')
    options = rmfield(options, 'right_component_fun');
end
if isfield(options, 'right_bands')
    options = rmfield(options, 'right_bands');
end
space.options = options;
space.solve_info = solve_info;
space.adaptive_info = struct('enabled', false, 'residual', NaN, ...
    'initial_rank', options.rank, 'final_rank', options.rank, ...
    'tolerance', NaN, 'reached_tolerance', false);
space.sampling_swapped_left_right = false;
space.ngrid = ngrid;
space.nleft = nleft;
space.nright = nright;
local_progress(options, 1, 'complete');
end

function weight = local_streamed_weight( ...
    left, provider, right_bands, block_size, options)
% Accumulate the same k-means weight without retaining all right bands.
ngrid = size(left{1}, 1);
effective_weight = options.weight;
if (~isfield(options, 'weight_was_set') || ~options.weight_was_set) && ...
        numel(left) > 1
    effective_weight = 'prod';
end
switch lower(effective_weight)
    case {'prod', 'multiply', 'power'}
        weight = zeros(ngrid, 1, 'like', left{1});
        for first = 1:block_size:numel(right_bands)
            bands = right_bands(first:min( ...
                numel(right_bands), first + block_size - 1));
            right = provider(bands);
            if ~iscell(right)
                right = {right};
            end
            for icomponent = 1:numel(left)
                left_i = left{icomponent};
                right_i = right{icomponent};
                for jcomponent = 1:numel(left)
                    left_j = left{jcomponent};
                    right_j = right{jcomponent};
                    weight = weight + sum(conj(left_i) .* left_j, 2) .* ...
                        sum(right_i .* conj(right_j), 2);
                end
            end
        end
        if strcmpi(effective_weight, 'power')
            weight = weight .^ (options.power / 2);
        end
    case 'add'
        weight = zeros(ngrid, 1, 'like', left{1});
        for icomponent = 1:numel(left)
            weight = weight + sum(abs(left{icomponent}).^2, 2);
        end
        for first = 1:block_size:numel(right_bands)
            bands = right_bands(first:min( ...
                numel(right_bands), first + block_size - 1));
            right = provider(bands);
            if ~iscell(right)
                right = {right};
            end
            for icomponent = 1:numel(left)
                weight = weight + sum(abs(right{icomponent}).^2, 2);
            end
        end
    otherwise
        error('ISDF:UnknownWeight', ...
            'Unknown ISDF weight "%s".', effective_weight);
end
weight = gather_if_gpu(real(weight(:)));
weight(~isfinite(weight) | weight < 0) = 0;
if all(weight == 0)
    weight = ones(size(weight));
end
end

function tf = local_materialize_product_mu(options, nleft, nright)
switch options.product_mu_storage
    case 'full'
        tf = true;
    case 'stream'
        tf = false;
    otherwise
        tf = options.rank * nleft * nright <= options.product_mu_max_elements;
end
end

function tf = local_stream_product_gram(options, ngrid)
switch options.product_gram_storage
    case 'full'
        tf = false;
    case 'stream'
        tf = true;
    otherwise
        tf = ngrid * options.rank > options.product_gram_max_elements;
end
end

function local_progress(options, fraction, stage)
if isfield(options, 'progress') && isa(options.progress, 'function_handle')
    options.progress(fraction, stage);
end
end

function local_print_rank(options, ngrid, nleft, nright)
isdf.report_rank(options, ngrid, nleft, nright);
end
