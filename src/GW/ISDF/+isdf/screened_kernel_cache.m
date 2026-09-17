function varargout = screened_kernel_cache(action, varargin)
%ISDF.SCREENED_KERNEL_CACHE Per-run store for projected reduced kernels.
%   ISDF.SCREENED_KERNEL_CACHE('reset') clears every stored kernel and the
%   hit/miss counters.
%
%   [VALUE, HIT] = ISDF.SCREENED_KERNEL_CACHE('get', KEY) returns a stored
%   kernel together with HIT = false when KEY is absent.
%
%   ISDF.SCREENED_KERNEL_CACHE('put', KEY, VALUE) stores VALUE when it fits
%   in the remaining byte budget.  Entries are never evicted: a sigma run
%   walks the same (k, q, target-space) keys in the same order for every
%   diagonal band, so keeping the entries that fit gives a stable hit
%   pattern, whereas an LRU policy would evict exactly the key needed next
%   and could thrash.  Values are shared and not copied, so callers must not
%   modify a kernel returned by 'get'.
%
%   INFO = ISDF.SCREENED_KERNEL_CACHE('stats') reports limit, bytes, stored,
%   hits, misses and skipped entries.
%
%   ISDF.SCREENED_KERNEL_CACHE('limit', BYTES) sets the byte budget used by
%   later 'put' calls.  A non-positive budget disables storage.

persistent initialized keys values sizes total_bytes limit hits misses skipped

if isempty(initialized)
    initialized = true;
    keys = {};
    values = {};
    sizes = zeros(1, 0);
    total_bytes = 0;
    limit = 0;
    hits = 0;
    misses = 0;
    skipped = 0;
end

switch lower(action)
    case 'reset'
        keys = {};
        values = {};
        sizes = zeros(1, 0);
        total_bytes = 0;
        hits = 0;
        misses = 0;
        skipped = 0;
    case 'limit'
        limit = local_limit(local_argument(varargin, 1, action, 'BYTES'));
    case 'get'
        key = local_argument(varargin, 1, action, 'KEY');
        index = find(strcmp(keys, key), 1);
        if isempty(index)
            misses = misses + 1;
            varargout = {[], false};
        else
            hits = hits + 1;
            varargout = {values{index}, true};
        end
    case 'put'
        key = local_argument(varargin, 1, action, 'KEY');
        value = local_argument(varargin, 2, action, 'VALUE');
        bytes = local_bytes(value);
        index = find(strcmp(keys, key), 1);
        if isempty(index)
            if bytes > limit - total_bytes
                skipped = skipped + 1;
                return;
            end
            keys{end + 1} = key;
            values{end + 1} = value;
            sizes(end + 1) = bytes;
            total_bytes = total_bytes + bytes;
        else
            if total_bytes - sizes(index) + bytes > limit
                skipped = skipped + 1;
                return;
            end
            values{index} = value;
            total_bytes = total_bytes - sizes(index) + bytes;
            sizes(index) = bytes;
        end
    case 'stats'
        varargout = {struct('limit', limit, 'bytes', total_bytes, ...
            'stored', numel(keys), 'hits', hits, 'misses', misses, ...
            'skipped', skipped)};
    otherwise
        error('ISDF:UnknownScreenedKernelCacheAction', ...
            'Unknown screened-kernel cache action "%s".', action);
end
end

function value = local_argument(args, index, action, name)
if numel(args) < index
    error('ISDF:ScreenedKernelCacheArgument', ...
        'ISDF.SCREENED_KERNEL_CACHE(''%s'', ...) requires the %s argument.', ...
        action, name);
end
value = args{index};
end

function bytes = local_limit(value)
if ~(isnumeric(value) && isscalar(value) && value >= 0)
    error('ISDF:ScreenedKernelCacheLimit', ...
        ['Screened-kernel cache limit must be a non-negative scalar ' ...
         '(Inf for unlimited).']);
end
bytes = double(value);
end

function bytes = local_bytes(value)
% Conservative element width: complex arrays count both parts, and GPU
% arrays are assumed complex because the stored projected kernels are.
if isa(value, 'gpuArray')
    complex_factor = 2;
else
    complex_factor = 1 + double(~isreal(value));
end
switch class(value)
    case 'single'
        width = 4 * complex_factor;
    otherwise
        width = 8 * complex_factor;
end
bytes = numel(value) * width;
end
