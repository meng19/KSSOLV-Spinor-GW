function [x, info] = stable_solve(a, b, options)
%STABLE_SOLVE Solve x*b = a with an ISDF-aware fallback.

if nargin < 3 || isempty(options)
    options = struct();
end
if ~isfield(options, 'rcond_tol') || isempty(options.rcond_tol)
    options.rcond_tol = 1e-12;
end
if ~isfield(options, 'warn_ill_conditioned') || isempty(options.warn_ill_conditioned)
    options.warn_ill_conditioned = false;
end

info = struct();
if size(b, 1) == size(b, 2)
    info.rcond = rcond(b);
else
    singular_values = svd(b, 'econ');
    if isempty(singular_values) || max(singular_values) == 0
        info.rcond = 0;
    else
        info.rcond = min(singular_values) / max(singular_values);
    end
end
info.used_pinv = false;
info.ill_conditioned = info.rcond < options.rcond_tol;

if info.ill_conditioned && options.warn_ill_conditioned
    warning('ISDF:IllConditionedInterpolation', ...
        ['ISDF interpolation matrix is ill-conditioned (rcond = %.3e). ' ...
         'Using direct solve to preserve reference-path behavior.'], info.rcond);
end
if isfield(options, 'ill_conditioned_solver') && ...
        strcmpi(options.ill_conditioned_solver, 'pinv') && info.ill_conditioned
    info.used_pinv = true;
    if options.warn_ill_conditioned
        warning('ISDF:IllConditionedInterpolation', ...
            ['ISDF interpolation matrix is ill-conditioned (rcond = %.3e). ' ...
             'Using pseudoinverse fallback.'], info.rcond);
    end
    x = a * pinv(b);
else
    if info.ill_conditioned && ~options.warn_ill_conditioned
        warning_state_near = warning('off', 'MATLAB:nearlySingularMatrix');
        warning_state_sing = warning('off', 'MATLAB:singularMatrix');
        cleanup = onCleanup(@() restore_warnings( ...
            warning_state_near, warning_state_sing)); %#ok<NASGU>
    end
    x = a / b;
end
end

function restore_warnings(warning_state_near, warning_state_sing)
warning(warning_state_near);
warning(warning_state_sing);
end
