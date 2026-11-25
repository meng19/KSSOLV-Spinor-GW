function [cache_hit, eps] = epsilon_checkpoint_restore( ...
    filename, request_eps, sys, options, syms)
%EPSILON_CHECKPOINT_RESTORE Reuse a compatible completed epsilon checkpoint.

cache_hit = false;
eps = [];
if exist(char(filename), 'file') ~= 2
    return;
end

[cached_eps, cached_sys, cached_options, cached_syms, cached_request] = ...
    epsilon_checkpoint_load(filename);
if isempty(cached_request)
    fprintf(['Ignoring epsilon checkpoint without cache metadata; ' ...
        'recalculating.\n']);
    return;
end

% request_eps contains eps.isdf, so direct and all ISDF configurations are
% separate cache entries even when they share a checkpoint filename.
request_match = isequaln(cached_request, request_eps);
sys_match = isequaln(cached_sys, sys);
options_match = isequaln(cached_options, options);
syms_match = isequaln(cached_syms, syms);

if request_match && sys_match && options_match && syms_match
    cache_hit = true;
    eps = cached_eps;
else
    fprintf('Ignoring incompatible epsilon checkpoint; recalculating.\n');
    local_report_difference('request_eps', cached_request, request_eps, request_match);
    local_report_difference('sys', cached_sys, sys, sys_match);
    local_report_difference('options', cached_options, options, options_match);
    local_report_difference('syms', cached_syms, syms, syms_match);
end
end

function local_report_difference(label, cached, current, matched)
%LOCAL_REPORT_DIFFERENCE Name the data that invalidated a cached result.
%   Diagnostic only: the cache decision above stays one isequaln comparison
%   per component, so this never turns a miss into a hit or vice versa.
if matched
    return;
end
if ~strcmp(class(cached), class(current))
    fprintf('  %s: class changed from %s to %s\n', ...
        label, class(cached), class(current));
    return;
end
cached_struct = local_as_struct(cached);
current_struct = local_as_struct(current);
if numel(cached_struct) ~= 1 || numel(current_struct) ~= 1
    fprintf('  %s: value differs\n', label);
    return;
end
names = union(fieldnames(cached_struct), fieldnames(current_struct));
changed = {};
for ii = 1:numel(names)
    name = names{ii};
    if ~isfield(cached_struct, name)
        changed{end+1} = [name ' (new)']; %#ok<AGROW>
    elseif ~isfield(current_struct, name)
        changed{end+1} = [name ' (removed)']; %#ok<AGROW>
    elseif ~isequaln(cached_struct.(name), current_struct.(name))
        changed = local_describe_change(name, ...
            cached_struct.(name), current_struct.(name), changed);
    end
end
if isempty(changed)
    fprintf(['  %s: object differs outside the public fields ' ...
        '(private property or class-specific equality)\n'], label);
else
    fprintf('  %s: %s\n', label, local_join(changed));
end
end

function changed = local_describe_change(name, cached, current, changed)
%LOCAL_DESCRIBE_CHANGE Describe one changed field, expanding one struct level.
if isstruct(cached) && isstruct(current) && isscalar(cached) && isscalar(current)
    sub_names = union(fieldnames(cached), fieldnames(current));
    for ii = 1:numel(sub_names)
        sub = sub_names{ii};
        sub_name = [name '.' sub];
        if ~isfield(cached, sub)
            changed{end+1} = [sub_name ' (new)']; %#ok<AGROW>
        elseif ~isfield(current, sub)
            changed{end+1} = [sub_name ' (removed)']; %#ok<AGROW>
        elseif ~isequaln(cached.(sub), current.(sub))
            changed{end+1} = local_value_change( ...
                sub_name, cached.(sub), current.(sub)); %#ok<AGROW>
        end
    end
else
    changed{end+1} = local_value_change(name, cached, current);
end
end

function text = local_value_change(name, cached, current)
%LOCAL_VALUE_CHANGE Append the values when they are scalar and cheap to print.
if local_is_scalar_number(cached) && local_is_scalar_number(current)
    text = sprintf('%s (%g -> %g)', name, double(cached), double(current));
elseif isrow(cached) && ischar(cached) && isrow(current) && ischar(current)
    text = sprintf('%s (%s -> %s)', name, cached, current);
else
    text = name;
end
end

function tf = local_is_scalar_number(value)
tf = (isnumeric(value) || islogical(value)) && isscalar(value);
end

function s = local_as_struct(value)
%LOCAL_AS_STRUCT Expose the public fields of a struct or a value-class object.
if isstruct(value)
    s = value;
elseif isobject(value)
    try
        s = struct(value);
    catch
        s = struct();
    end
else
    s = struct();
end
end

function text = local_join(items)
%LOCAL_JOIN Join changed names, truncating long lists.
limit = 12;
if numel(items) > limit
    items = [items(1:limit), {sprintf('... (%d more)', numel(items) - limit)}];
end
text = strjoin(items, ', ');
end
