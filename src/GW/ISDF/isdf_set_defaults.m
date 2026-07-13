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
end
