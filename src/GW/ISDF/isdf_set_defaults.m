function options = isdf_set_defaults(options, nphi, npsi, ngrid)
%ISDF_SET_DEFAULTS Fill a small, local ISDF option structure.

if nargin < 1 || isempty(options)
    options = struct();
end

max_rank = min([ngrid, nphi * npsi]);
if ~isfield(options, 'rank') || isempty(options.rank)
    options.rank = max_rank;
end
options.rank = min(max(1, ceil(options.rank)), max_rank);

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
