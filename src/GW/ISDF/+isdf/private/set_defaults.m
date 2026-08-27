function options = set_defaults(options, nleft, nright, ngrid)
%SET_DEFAULTS Fill ISDF product-space option defaults.

max_rank = min([ngrid, nleft * nright]);
rank_was_set = isfield(options, 'rank') && ~isempty(options.rank);
if ~isfield(options, 'rank_ratio') || isempty(options.rank_ratio)
    options.rank_ratio = 1;
end
if ~isfield(options, 'warn_rank_selection') || isempty(options.warn_rank_selection)
    options.warn_rank_selection = true;
end
recommended_rank = ceil(sqrt(nleft * nright) * options.rank_ratio);
recommended_rank = min(max(1, recommended_rank), max_rank);
options.max_rank = max_rank;
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
if ~isfield(options, 'rcond_tol') || isempty(options.rcond_tol)
    options.rcond_tol = 1e-12;
end
if ~isfield(options, 'warn_ill_conditioned') || isempty(options.warn_ill_conditioned)
    options.warn_ill_conditioned = false;
end
if ~isfield(options, 'weight') || isempty(options.weight)
    options.weight = 'add';
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
if ~isfield(options, 'random_oversampling') || isempty(options.random_oversampling)
    options.random_oversampling = 1.2;
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
