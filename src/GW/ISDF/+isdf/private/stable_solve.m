function [x, info] = stable_solve(a, b, options)
%STABLE_SOLVE Solve x*b = a with an ISDF-aware fallback.

info = struct();
info.used_pinv = false;
compute_rcond = ~isfield(options, 'compute_rcond') || ...
    options.compute_rcond;
needs_condition = compute_rcond || options.warn_ill_conditioned || ...
    (isfield(options, 'ill_conditioned_solver') && ...
    strcmpi(options.ill_conditioned_solver, 'pinv'));
if ~needs_condition
    info.rcond = NaN;
    info.ill_conditioned = false;
    x = a / b;
    return;
end
b_condition = gather_if_gpu(b);
if size(b, 1) == size(b, 2)
    info.rcond = rcond(b_condition);
else
    singular_values = svd(b_condition, 'econ');
    if isempty(singular_values) || max(singular_values) == 0
        info.rcond = 0;
    else
        info.rcond = min(singular_values) / max(singular_values);
    end
end
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
