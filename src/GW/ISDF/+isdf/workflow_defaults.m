function options = workflow_defaults(options, freq_dep)
%WORKFLOW_DEFAULTS Fill defaults shared by epsilon and sigma ISDF workflows.
%
% Product-space and numerical defaults belong to private/set_defaults and
% are applied when isdf.build_space or isdf.polarizability is entered.

if nargin < 1 || isempty(options) || ~isstruct(options)
    options = struct();
end
if nargin < 2 || isempty(freq_dep)
    freq_dep = 0;
end

if ~isfield(options, 'enable') || isempty(options.enable)
    options.enable = false;
end
if ~isfield(options, 'algorithm') || isempty(options.algorithm)
    if freq_dep == 0
        options.algorithm = 'reduced_basis';
    else
        options.algorithm = 'matrix_elements';
    end
end
if ~isfield(options, 'validate_hf_exchange') || ...
        isempty(options.validate_hf_exchange)
    % Direct FFT validation repeats selected product transforms.
    options.validate_hf_exchange = false;
end
if ~isfield(options, 'validate_hf_exchange_max_pairs') || ...
        isempty(options.validate_hf_exchange_max_pairs)
    options.validate_hf_exchange_max_pairs = Inf;
end
if ~isfield(options, 'reduced_solver') || isempty(options.reduced_solver)
    options.reduced_solver = 'cauchy';
end
if ~isfield(options, 'cauchy_froErr') || isempty(options.cauchy_froErr)
    options.cauchy_froErr = 1e-8;
end
if ~isfield(options, 'cauchy_MaxIter') || isempty(options.cauchy_MaxIter)
    options.cauchy_MaxIter = 12;
end
end
