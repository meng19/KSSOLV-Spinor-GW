function varargout = sigma_isdf_component_cache(action, varargin)
%SIGMA_ISDF_COMPONENT_CACHE Bounded LRU cache for sigma ISDF spaces/components.

persistent keys values sizes limit total_bytes;
if isempty(limit)
    limit = Inf;
    total_bytes = 0;
end
switch lower(action)
    case 'reset'
        keys = {};
        values = {};
        sizes = [];
        total_bytes = 0;
    case 'limit'
        limit = varargin{1};
        while total_bytes > limit && ~isempty(keys)
            local_remove(1);
        end
    case 'get'
        key = varargin{1};
        index = find(strcmp(keys, key), 1);
        if isempty(index)
            varargout = {[], false};
        else
            value = values{index};
            byte = sizes(index);
            local_remove(index);
            keys{end + 1} = key;
            values{end + 1} = value;
            sizes(end + 1) = byte;
            varargout = {value, true};
        end
    case 'put'
        key = varargin{1};
        value = varargin{2};
        if numel(varargin) >= 3
            bytes = varargin{3};
        else
            bytes = local_estimate_bytes(value);
        end
        index = find(strcmp(keys, key), 1);
        if ~isempty(index)
            local_remove(index);
        end
        if isfinite(limit) && bytes > limit
            return;
        end
        while isfinite(limit) && total_bytes + bytes > limit && ~isempty(keys)
            local_remove(1);
        end
        keys{end + 1} = key;
        values{end + 1} = value;
        sizes(end + 1) = bytes;
        total_bytes = total_bytes + bytes;
    otherwise
        error('Sigma:UnknownISDFCacheAction', ...
            'Unknown ISDF component cache action "%s".', action);
end

    function local_remove(index)
        total_bytes = total_bytes - sizes(index);
        keys(index) = [];
        values(index) = [];
        sizes(index) = [];
    end
end

function bytes = local_estimate_bytes(value)
% Conservative estimate; shared arrays may be counted more than once.
if isnumeric(value) || islogical(value)
    bytes = 16 * numel(value);
elseif isstruct(value)
    bytes = 0;
    fields = fieldnames(value);
    for ielement = 1:numel(value)
        for ifield = 1:numel(fields)
            bytes = bytes + local_estimate_bytes( ...
                value(ielement).(fields{ifield}));
        end
    end
elseif iscell(value)
    bytes = 0;
    for ielement = 1:numel(value)
        bytes = bytes + local_estimate_bytes(value{ielement});
    end
else
    bytes = 0;
end
end
