function varargout = sigma_isdf_component_cache(action, varargin)
%SIGMA_ISDF_COMPONENT_CACHE Per-sigma-run cache for right real components.

persistent keys values;
switch lower(action)
    case 'reset'
        keys = {};
        values = {};
    case 'get'
        key = varargin{1};
        index = find(strcmp(keys, key), 1);
        if isempty(index)
            varargout = {[], false};
        else
            varargout = {values{index}, true};
        end
    case 'put'
        key = varargin{1};
        value = varargin{2};
        keys{end + 1} = key;
        values{end + 1} = value;
    otherwise
        error('Sigma:UnknownISDFCacheAction', ...
            'Unknown ISDF component cache action "%s".', action);
end
end
