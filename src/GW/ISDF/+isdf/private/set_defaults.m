function options = set_defaults(options, nleft, nright, ngrid, fftgrid)
%SET_DEFAULTS Fill ISDF product-space option defaults.

if nargin < 1 || isempty(options)
    options = struct();
end
if nargin < 2
    options = local_polarizability_defaults(options);
    return;
end
if nargin < 5
    fftgrid = [];
end

max_rank = min([ngrid, nleft * nright]);
rank_was_set = isfield(options, 'rank') && ~isempty(options.rank);
if ~isfield(options, 'rank_ratio') || isempty(options.rank_ratio)
    options.rank_ratio = 1;
end
if ~isfield(options, 'warn_rank_selection') || isempty(options.warn_rank_selection)
    options.warn_rank_selection = true;
end
unclamped_recommended_rank = ceil(sqrt(nleft * nright) * options.rank_ratio);
recommended_rank = min(max(1, unclamped_recommended_rank), max_rank);
options.max_rank = max_rank;
options.unclamped_recommended_rank = unclamped_recommended_rank;
options.recommended_rank = recommended_rank;
if ~isfield(options, 'rank') || isempty(options.rank)
    options.rank = recommended_rank;
    options.rank_source = 'auto';
else
    requested_rank = options.rank;
    options.rank_source = 'user';
end
options.rank = min(max(1, ceil(options.rank)), max_rank);
if rank_was_set && options.warn_rank_selection
    if options.rank ~= ceil(requested_rank)
        warn_once('ISDF:RankClamped', ...
            sprintf('%.16g:%d:%d', requested_rank, options.rank, max_rank), ...
            ['ISDF rank %.3g was clamped to %d. Allowed range is ' ...
             '1 <= rank <= %d.'], requested_rank, options.rank, max_rank);
    elseif options.rank > recommended_rank
        warn_once('ISDF:RankAboveRecommended', ...
            sprintf('%d:%d:%d:%d:%.16g', options.rank, recommended_rank, ...
            nleft, nright, options.rank_ratio), ...
            ['ISDF rank %d exceeds the recommended default %d ' ...
             '(ceil(sqrt(%d*%d)*rank_ratio), rank_ratio = %.3g). ' ...
             'Large ranks can make the interpolation solve ill-conditioned; ' ...
             'keep this value only for convergence tests or full-rank ' ...
             'validation.'], options.rank, recommended_rank, nleft, nright, ...
             options.rank_ratio);
    elseif options.rank < recommended_rank
        warn_once('ISDF:RankBelowRecommended', ...
            sprintf('%d:%d:%d:%d:%.16g', options.rank, recommended_rank, ...
            nleft, nright, options.rank_ratio), ...
            ['ISDF rank %d is below the recommended default %d ' ...
             '(ceil(sqrt(%d*%d)*rank_ratio), rank_ratio = %.3g). ' ...
             'Small ranks can under-resolve the product space and increase ' ...
             'matrix-element errors; keep this value only for convergence ' ...
             'tests or deliberately compressed runs.'], options.rank, ...
             recommended_rank, nleft, nright, options.rank_ratio);
    end
end
if ~isfield(options, 'sample_method') || isempty(options.sample_method)
    options.sample_method = 'qrcp';
end
if ~isfield(options, 'seed') || isempty(options.seed)
    options.seed = 0;
end
if ~isfield(options, 'swap_left_right') || isempty(options.swap_left_right)
    % Select points from conj(right).*left, then reconstruct the original
    % conj(left).*right space.  This is a validation control only.
    options.swap_left_right = false;
end
if ~(isscalar(options.swap_left_right) && ...
        (islogical(options.swap_left_right) || isnumeric(options.swap_left_right)))
    error('ISDF:SwapLeftRight', 'swap_left_right must be a logical scalar.');
end
options.swap_left_right = logical(options.swap_left_right);
for product_type = {'vc', 'vn', 'nn'}
    suffix = product_type{1};
    for base_name = {'rank', 'rank_ratio', 'swap_left_right'}
        field = [base_name{1} '_' suffix];
        if ~isfield(options, field)
            options.(field) = [];
        end
    end
end
if ~isfield(options, 'rcond_tol') || isempty(options.rcond_tol)
    options.rcond_tol = 1e-12;
end
if ~isfield(options, 'warn_ill_conditioned') || isempty(options.warn_ill_conditioned)
    options.warn_ill_conditioned = false;
end
if ~isfield(options, 'ill_conditioned_solver') || ...
        isempty(options.ill_conditioned_solver)
    options.ill_conditioned_solver = 'direct';
end
options.ill_conditioned_solver = lower(char(options.ill_conditioned_solver));
if ~any(strcmp(options.ill_conditioned_solver, {'direct', 'pinv'}))
    error('ISDF:IllConditionedSolver', ...
        'ill_conditioned_solver must be ''direct'' or ''pinv''.');
end
if ~isfield(options, 'weight') || isempty(options.weight)
    options.weight = 'add';
end
if ~isfield(options, 'weight_was_set') || isempty(options.weight_was_set)
    options.weight_was_set = false;
end
if ~isfield(options, 'points') || isempty(options.points)
    options.points = [];
end
if ~isfield(options, 'progress') || isempty(options.progress)
    options.progress = [];
end
if ~isfield(options, 'print_rank') || isempty(options.print_rank)
    options.print_rank = true;
end
if ~isfield(options, 'product_type') || isempty(options.product_type)
    options.product_type = '';
end
if ~isfield(options, 'fftgrid') || isempty(options.fftgrid)
    options.fftgrid = fftgrid;
end
if ~isfield(options, 'validate_hf_exchange_max_pairs') || ...
        isempty(options.validate_hf_exchange_max_pairs)
    options.validate_hf_exchange_max_pairs = Inf;
end
if ~isfield(options, 'validate_hf_exchange_label') || ...
        isempty(options.validate_hf_exchange_label)
    options.validate_hf_exchange_label = 'HF exchange';
end
if ~isfield(options, 'power') || isempty(options.power)
    options.power = 1;
end
if ~isfield(options, 'init') || isempty(options.init)
    options.init = 'random';
end
if ~isfield(options, 'kmeans_max_iter') || isempty(options.kmeans_max_iter)
    options.kmeans_max_iter = 100;
end
if ~isfield(options, 'kmeans_distance_block_elements') || ...
        isempty(options.kmeans_distance_block_elements)
    % Default to the full distance matrix for maximum BLAS throughput.
    % Memory-limited calculations must set a finite block-element limit.
    options.kmeans_distance_block_elements = Inf;
end
if ~(isscalar(options.kmeans_distance_block_elements) && ...
        options.kmeans_distance_block_elements >= 1)
    error('ISDF:KmeansDistanceBlockElements', ...
        'kmeans_distance_block_elements must be a positive scalar.');
end
if ~isfield(options, 'random_oversampling') || isempty(options.random_oversampling)
    options.random_oversampling = 1.2;
end
if ~isfield(options, 'sample_precision') || isempty(options.sample_precision)
    options.sample_precision = 'double';
end
options.sample_precision = lower(char(options.sample_precision));
if ~any(strcmp(options.sample_precision, {'double', 'single'}))
    error('ISDF:SamplePrecision', ...
        'sample_precision must be ''double'' or ''single'' (got ''%s'').', ...
        options.sample_precision);
end
if ~isfield(options, 'projection_block_elements') || ...
        isempty(options.projection_block_elements)
    % Do not split projected product calculations unless requested.
    options.projection_block_elements = Inf;
end
if ~(isscalar(options.projection_block_elements) && ...
        options.projection_block_elements >= 1)
    error('ISDF:ProjectionBlockElements', ...
        'projection_block_elements must be a positive scalar.');
end
if ~isfield(options, 'product_mu_storage') || isempty(options.product_mu_storage)
    % The dense coefficient matrix is the fastest representation.
    options.product_mu_storage = 'full';
end
options.product_mu_storage = lower(char(options.product_mu_storage));
if ~any(strcmp(options.product_mu_storage, {'auto', 'full', 'stream'}))
    error('ISDF:ProductMuStorage', ...
        'product_mu_storage must be ''auto'', ''full'', or ''stream''.');
end
if ~isfield(options, 'product_mu_max_elements') || ...
        isempty(options.product_mu_max_elements)
    % Retain the dense, faster representation through roughly 1.5 GiB for
    % complex double data.  Larger NN spaces switch to streamed coefficients.
    options.product_mu_max_elements = 1e8;
end
if ~(isscalar(options.product_mu_max_elements) && ...
        options.product_mu_max_elements >= 1)
    error('ISDF:ProductMuMaxElements', ...
        'product_mu_max_elements must be a positive scalar.');
end
if ~isfield(options, 'product_gram_storage') || ...
        isempty(options.product_gram_storage)
    % The full C1 path avoids block-loop and repeated-FFT overhead.
    options.product_gram_storage = 'full';
end
options.product_gram_storage = lower(char(options.product_gram_storage));
if ~any(strcmp(options.product_gram_storage, {'auto', 'full', 'stream'}))
    error('ISDF:ProductGramStorage', ...
        'product_gram_storage must be ''auto'', ''full'', or ''stream''.');
end
if ~isfield(options, 'product_gram_max_elements') || ...
        isempty(options.product_gram_max_elements)
    % Keep the full C1 path by default whenever it is reasonably sized;
    % streaming remains the automatic safeguard for large grids/ranks.
    options.product_gram_max_elements = 1e8;
end
if ~(isscalar(options.product_gram_max_elements) && ...
        options.product_gram_max_elements >= 1)
    error('ISDF:ProductGramMaxElements', ...
        'product_gram_max_elements must be a positive scalar.');
end
if ~isfield(options, 'product_gram_block_elements') || ...
        isempty(options.product_gram_block_elements)
    if strcmp(options.product_gram_storage, 'stream')
        % A caller explicitly requesting streaming expects a bounded block.
        options.product_gram_block_elements = 2e6;
    else
        options.product_gram_block_elements = Inf;
    end
end
if ~(isscalar(options.product_gram_block_elements) && ...
        options.product_gram_block_elements >= 1)
    error('ISDF:ProductGramBlockElements', ...
        'product_gram_block_elements must be a positive scalar.');
end
if ~isfield(options, 'interpolation_solver') || isempty(options.interpolation_solver)
    options.interpolation_solver = 'direct';
end
options.interpolation_solver = lower(char(options.interpolation_solver));
if ~any(strcmp(options.interpolation_solver, {'direct', 'svd_whiten'}))
    error('ISDF:InterpolationSolver', ...
        'interpolation_solver must be ''direct'' or ''svd_whiten''.');
end
if ~isfield(options, 'svd_cutoff') || isempty(options.svd_cutoff)
    options.svd_cutoff = 0;
end
if ~isfield(options, 'svd_ratio') || isempty(options.svd_ratio)
    options.svd_ratio = 0.5;
end
if options.svd_cutoff < 0 || options.svd_cutoff >= 1 || ...
        options.svd_ratio < 0 || options.svd_ratio > 1
    error('ISDF:SVDOptions', ...
        'svd_cutoff must be in [0,1), and svd_ratio must be in [0,1].');
end
if ~isfield(options, 'adaptive_rank_enable') || isempty(options.adaptive_rank_enable)
    options.adaptive_rank_enable = false;
end
if ~isfield(options, 'adaptive_rank_tol') || isempty(options.adaptive_rank_tol)
    options.adaptive_rank_tol = 1e-6;
end
if ~isfield(options, 'adaptive_rank_step') || isempty(options.adaptive_rank_step)
    options.adaptive_rank_step = max(1, ceil(sqrt(nleft * nright)));
end
if ~isfield(options, 'adaptive_rank_max') || isempty(options.adaptive_rank_max)
    options.adaptive_rank_max = max_rank;
end
if ~isfield(options, 'adaptive_validation_rank') || ...
        isempty(options.adaptive_validation_rank)
    options.adaptive_validation_rank = max(16, options.adaptive_rank_step);
end
options.adaptive_rank_max = min(max(1, ceil(options.adaptive_rank_max)), max_rank);
options.adaptive_rank_step = max(1, ceil(options.adaptive_rank_step));
options.adaptive_validation_rank = max(1, ...
    ceil(options.adaptive_validation_rank));
end

function options = local_polarizability_defaults(options)
% Keep reduced-polarizability defaults in the same central ISDF entry point.
if ~isfield(options, 'method') || isempty(options.method)
    options.method = 'cauchy';
end
if ~isfield(options, 'freq') || isempty(options.freq)
    options.freq = 0;
end
if ~isfield(options, 'progress') || isempty(options.progress)
    options.progress = [];
end
if ~isfield(options, 'froErr') || isempty(options.froErr)
    options.froErr = 1e-8;
end
if ~isfield(options, 'MaxIter') || isempty(options.MaxIter)
    options.MaxIter = 12;
end
if ~isfield(options, 'cauchy_coefficient_storage') || ...
        isempty(options.cauchy_coefficient_storage)
    options.cauchy_coefficient_storage = 'full';
end
options.cauchy_coefficient_storage = lower( ...
    char(options.cauchy_coefficient_storage));
if ~any(strcmp(options.cauchy_coefficient_storage, {'auto', 'full', 'stream'}))
    error('ISDF:CauchyCoefficientStorage', ...
        'cauchy_coefficient_storage must be ''auto'', ''full'', or ''stream''.');
end
if ~isfield(options, 'cauchy_coefficient_cache_bytes') || ...
        isempty(options.cauchy_coefficient_cache_bytes)
    options.cauchy_coefficient_cache_bytes = 4 * 1024^3;
end
if ~isfield(options, 'cauchy_coefficient_block_bytes') || ...
        isempty(options.cauchy_coefficient_block_bytes)
    if strcmp(options.cauchy_coefficient_storage, 'full')
        options.cauchy_coefficient_block_bytes = Inf;
    else
        options.cauchy_coefficient_block_bytes = 512 * 1024^2;
    end
end
if ~(isscalar(options.cauchy_coefficient_cache_bytes) && ...
        options.cauchy_coefficient_cache_bytes >= 1 && ...
        isscalar(options.cauchy_coefficient_block_bytes) && ...
        options.cauchy_coefficient_block_bytes >= 1)
    error('ISDF:CauchyCoefficientBytes', ...
        'Cauchy coefficient byte limits must be positive scalars.');
end
end

function warn_once(identifier, key, message, varargin)
persistent issued_keys;
if isempty(issued_keys)
    issued_keys = {};
end
full_key = [identifier ':' key];
if any(strcmp(issued_keys, full_key))
    return;
end
issued_keys{end + 1} = full_key;
warning(identifier, message, varargin{:});
end
