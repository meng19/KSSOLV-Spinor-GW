function [ind_mu, products, adaptive_state] = sample_points(left, right, options)
%SAMPLE_POINTS Select ISDF interpolation points for component products.

if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
rng(options.seed, 'twister');
products = [];
adaptive_state = struct('enabled', false);
adaptive_methods = {'qrcp', 'qrcp_randomized', 'randomized_qrcp', 'default'};
if options.adaptive_rank_enable && ...
        ~any(strcmpi(options.sample_method, adaptive_methods))
    error('ISDF:AdaptiveRankRequiresQRCP', ...
        ['Adaptive rank selection requires QRCP sampling. Use ''qrcp'' ', ...
         'for an exact residual or ''qrcp_randomized'' for a projected ', ...
         'validation residual.']);
end

switch lower(options.sample_method)
    case 'qrcp'
        products = component_products(left, right, [], []);
        if options.adaptive_rank_enable
            [ind_mu, adaptive_state] = adaptive_qrcp_sample(products, options);
        else
            ind_mu = qrcp_sample(products, options.rank);
        end
    case {'qrcp_randomized', 'randomized_qrcp', 'default'}
        if options.adaptive_rank_enable
            [ind_mu, adaptive_state] = adaptive_randomized_qrcp_sample( ...
                left, right, options);
        elseif numel(left) == 1
            ind_mu = scalar_randomized_sample( ...
                left{1}, right{1}, options);
        else
            npairs = size(left{1}, 2) * size(right{1}, 2);
            projection_rank = max(options.rank, ...
                ceil(options.random_oversampling * options.rank));
            projection_rank = min(projection_rank, npairs);
            projection = randn_like(npairs, projection_rank, left{1});
            has_complex = any(cellfun(@(x) ~isreal(x), left)) || ...
                any(cellfun(@(x) ~isreal(x), right));
            if has_complex
                projection = projection + 1i * ...
                    randn_like(npairs, projection_rank, left{1});
            end
            projection = sample_cast(projection, options);
            compressed_products = component_products( ...
                left, right, [], projection, options.sample_precision);
            ind_mu = qrcp_sample(compressed_products, options.rank);
        end
    case 'kmeans'
        if numel(left) == 1
            weight = scalar_weight(left{1}, right{1}, options);
        else
            weight = component_weight(left, right, options);
        end
        ind_mu = kmeans_sample(weight, options);
    otherwise
        error('ISDF:UnknownSampleMethod', ...
            ['Unknown ISDF sample_method "%s". Supported methods: qrcp, ' ...
             'qrcp_randomized, kmeans, default.'], options.sample_method);
end

function [ind_mu, state] = adaptive_randomized_qrcp_sample( ...
    left, right, options)
%ADAPTIVE_RANDOMIZED_QRCP_SAMPLE QRCP rank selection with held-out sketches.

npairs = size(left{1}, 2) * size(right{1}, 2);
max_rank = min(max(options.rank, options.adaptive_rank_max), npairs);
train_rank = min(npairs, max(max_rank, ...
    ceil(options.random_oversampling * max_rank)));
validation_rank = min(npairs, options.adaptive_validation_rank);
train_projection = local_projection( ...
    npairs, train_rank, left, right, options);
validation_projection = local_projection( ...
    npairs, validation_rank, left, right, options);
train_products = component_products( ...
    left, right, [], train_projection, options.sample_precision);
validation_products = component_products( ...
    left, right, [], validation_projection, options.sample_precision);
all_pivots = qrcp_sample(train_products, max_rank);
validation_norm = gather_if_gpu(norm(validation_products, 'fro'));

initial_rank = options.rank;
if validation_norm == 0
    ind_mu = all_pivots(1:initial_rank);
    state = local_adaptive_state(initial_rank, initial_rank, 0, options, ...
        false, [], struct(), 'randomized_validation');
    return;
end

rank_mu = initial_rank;
residual = inf;
zeta = [];
solve_info = struct();
while true
    ind_mu = all_pivots(1:rank_mu);
    [c1, c2] = product_gram(left, right, ind_mu);
    [zeta, solve_info] = stable_solve(c1, c2, options);
    validation_mu = validation_products(ind_mu, :);
    residual = gather_if_gpu(norm( ...
        validation_products - zeta * validation_mu, 'fro') / validation_norm);
    if residual <= options.adaptive_rank_tol || rank_mu >= max_rank
        break;
    end
    rank_mu = min(max_rank, rank_mu + options.adaptive_rank_step);
end
state = local_adaptive_state(initial_rank, rank_mu, residual, options, ...
    true, zeta, solve_info, 'randomized_validation');
end

function projection = local_projection(npairs, ncolumns, left, right, options)
projection = randn_like(npairs, ncolumns, left{1});
has_complex = any(cellfun(@(x) ~isreal(x), left)) || ...
    any(cellfun(@(x) ~isreal(x), right));
if has_complex
    projection = projection + 1i * randn_like(npairs, ncolumns, left{1});
end
projection = sample_cast(projection, options);
end
end

function [ind_mu, state] = adaptive_qrcp_sample(products, options)
%ADAPTIVE_QRCP_SAMPLE Increase rank along QRCP pivots until residual converges.

initial_rank = options.rank;
max_rank = max(initial_rank, options.adaptive_rank_max);
max_rank = min(max_rank, min(size(products)));
all_pivots = qrcp_sample(products, max_rank);
if gather_if_gpu(norm(products, 'fro')) == 0
    ind_mu = all_pivots(1:initial_rank);
    state = local_adaptive_state(initial_rank, initial_rank, 0, options, ...
        false, [], struct(), 'exact_residual');
    return;
end

rank_mu = initial_rank;
residual = inf;
zeta = [];
solve_info = struct();
while true
    ind_mu = all_pivots(1:rank_mu);
    [residual, zeta, solve_info] = qrcp_interpolation_residual( ...
        products, ind_mu, options);
    if residual <= options.adaptive_rank_tol || rank_mu >= max_rank
        break;
    end
    rank_mu = min(max_rank, rank_mu + options.adaptive_rank_step);
end
state = local_adaptive_state(initial_rank, rank_mu, residual, options, ...
    true, zeta, solve_info, 'exact_residual');
end

function state = local_adaptive_state(initial_rank, final_rank, residual, ...
    options, has_zeta, zeta, solve_info, residual_kind)
state = struct('enabled', true, 'residual', residual, ...
    'initial_rank', initial_rank, 'final_rank', final_rank, ...
    'tolerance', options.adaptive_rank_tol, ...
    'reached_tolerance', residual <= options.adaptive_rank_tol, ...
    'residual_kind', residual_kind, 'has_zeta', has_zeta, ...
    'zeta', zeta, 'solve_info', solve_info);
end
