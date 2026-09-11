function method = gw_resolve_method(isdf_options, workflow)
%GW_RESOLVE_METHOD Resolve one canonical GW workflow method.

workflow = lower(char(workflow));
if ~any(strcmp(workflow, {'epsilon', 'sigma'}))
    error('ISDF:UnknownWorkflow', 'Unknown GW workflow "%s".', workflow);
end

if nargin < 1 || isempty(isdf_options) || ~isstruct(isdf_options) || ...
        ~isfield(isdf_options, 'enable') || ~isdf_options.enable
    method = 'direct';
    return;
end
if ~isfield(isdf_options, 'algorithm') || isempty(isdf_options.algorithm)
    error_id = sprintf('ISDF:Unknown%sAlgorithm', ...
        [upper(workflow(1)), workflow(2:end)]);
    error(error_id, '%s ISDF algorithm is missing.', workflow);
end

method = lower(char(isdf_options.algorithm));
if ~any(strcmp(method, {'matrix_elements', 'reduced_basis'}))
    error_id = sprintf('ISDF:Unknown%sAlgorithm', ...
        [upper(workflow(1)), workflow(2:end)]);
    error(error_id, ...
        ['Unknown %s ISDF algorithm "%s". Supported algorithms: ' ...
        'reduced_basis, matrix_elements.'], workflow, method);
end
end
