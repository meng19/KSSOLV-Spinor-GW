function [result, info] = comega_cstar(left, right, ev_occ, ev_unocc, options)
%COMEGA_CSTAR Dispatch reduced polarizability coefficient solvers.

if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
if ~isfield(options, 'method') || isempty(options.method)
    options.method = 'cauchy';
end
if ~isfield(options, 'freq') || isempty(options.freq)
    options.freq = 0;
end
ev_occ = ev_occ(:);
ev_unocc = ev_unocc(:);
freq = options.freq(:).';

switch lower(options.method)
    case 'direct'
        result = comega_direct(left, right, ev_occ, ev_unocc, freq);
        info = struct('method', 'direct', 'iterations', 0, ...
            'relative_error', 0);
    case 'cauchy'
        if numel(freq) ~= 1 || freq ~= 0
            error('ISDF:CauchyFullFrequencyUnsupported', ...
                ['Cauchy reduced polarizability currently supports ' ...
                 'static frequency only. Use reduced_solver = direct.']);
        end
        if ~isfield(options, 'froErr')
            options.froErr = 1e-8;
        end
        if ~isfield(options, 'MaxIter')
            options.MaxIter = 12;
        end
        [result, relative_error, iterations] = comega_cauchy( ...
            left, right, ev_occ, ev_unocc, options);
        info = struct('method', 'cauchy', ...
            'iterations', iterations, 'relative_error', relative_error);
    otherwise
        error('ISDF:UnknownCOmegaMethod', ...
            'Unknown COmegaCstar method "%s".', options.method);
end
end
