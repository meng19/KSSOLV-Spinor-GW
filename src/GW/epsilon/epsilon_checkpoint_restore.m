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
if isequaln(cached_request, request_eps) && isequaln(cached_sys, sys) && ...
        isequaln(cached_options, options) && isequaln(cached_syms, syms)
    cache_hit = true;
    eps = cached_eps;
else
    fprintf('Ignoring incompatible epsilon checkpoint; recalculating.\n');
end
end
