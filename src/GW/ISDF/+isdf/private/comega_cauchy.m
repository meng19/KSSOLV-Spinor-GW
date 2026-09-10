function [result, relative_error, iteration, fallback_direct] = comega_cauchy( ...
    left, right, ev_occ, ev_unocc, options)
%COMEGA_CAUCHY Cauchy reduced polarizability for frequency pages.

ncomponents = numel(left);
if ncomponents ~= numel(right)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end
if ~isfield(options, 'freq') || isempty(options.freq)
    options.freq = 0;
end
freq = options.freq(:).';
nmu = size(left{1}, 1);
nfreq = numel(freq);
result = complex(zeros(nmu, nmu, nfreq, 'like', left{1}));
relative_errors = zeros(1, nfreq);
iterations = zeros(1, nfreq);
fallback_direct = false(1, nfreq);
for ifreq = 1:nfreq
    omega = freq(ifreq);
    if omega == 0
        [result(:, :, ifreq), relative_errors(ifreq), iterations(ifreq)] = ...
            local_cauchy_resolvent(left, right, ev_occ, ev_unocc, options);
        continue;
    end

    [plus_page, plus_error, plus_iter, plus_ok] = ...
        local_cauchy_resolvent( ...
        left, right, ev_occ, ev_unocc - omega, options);
    [minus_page, minus_error, minus_iter, minus_ok] = ...
        local_cauchy_resolvent( ...
        left, right, ev_occ, ev_unocc + omega, options);
    if plus_ok && minus_ok
        result(:, :, ifreq) = 0.5 * (plus_page + minus_page);
        relative_errors(ifreq) = max(plus_error, minus_error);
        iterations(ifreq) = max(plus_iter, minus_iter);
    else
        direct_page = comega_direct(left, right, ev_occ, ev_unocc, omega);
        result(:, :, ifreq) = direct_page(:, :, 1);
        relative_errors(ifreq) = 0;
        iterations(ifreq) = 0;
        fallback_direct(ifreq) = true;
    end
end
relative_error = max(relative_errors);
iteration = max(iterations);
end

function [result, relative_error, iteration, ok] = local_cauchy_resolvent( ...
    left, right, ev_occ, ev_unocc, options)
ok = true;
[center, radius, ok] = local_contour(ev_occ, ev_unocc);
nmu = size(left{1}, 1);
if ~ok
    result = complex(zeros(nmu, nmu, 'like', left{1}));
    relative_error = inf;
    iteration = 0;
    return;
end
ev_occ_work = ev_occ;
ev_unocc_work = ev_unocc;
if isa(left{1}, 'gpuArray')
    ev_occ_work = gpuArray(ev_occ_work);
    ev_unocc_work = gpuArray(ev_unocc_work);
end
previous = [];
relative_error = inf;
for iteration = 1:options.MaxIter
    npoints = 2^(iteration + 3);
    if isempty(previous)
        result = complex(zeros(nmu, nmu, 'like', left{1}));
        point_indices = 0:npoints-1;
    else
        % The nodes for npoints/2 are the even nodes for npoints.  Reuse
        % their completed trapezoidal sum and evaluate only the new odd
        % nodes, rather than recomputing all Cauchy products.
        result = 0.5 * previous;
        point_indices = 1:2:npoints-1;
    end
    for ipoint = point_indices
        theta = 2 * pi * ipoint / npoints;
        exp_theta = exp(1i * theta);
        z = center + radius * exp_theta;
        occ_weight = 1 ./ (z - ev_occ_work);
        unocc_weight = 1 ./ (z - ev_unocc_work);
        result = result + local_weighted_products( ...
            left, right, occ_weight, unocc_weight) * ...
            (radius * exp_theta / npoints);
    end
    if ~isempty(previous)
        relative_error = gather_if_gpu( ...
            norm(result - previous, 'fro') / max(1, norm(result, 'fro')));
        if relative_error <= options.froErr
            return;
        end
    end
    previous = result;
end
end

function value = local_weighted_products(left, right, occ_weight, unocc_weight)
% Exploit the separable pair weight w(v,c)=w_v(v)*w_c(c).  Forming the
% full Nmu-by-(Nv*Nc) product matrix here is both the dominant allocation
% and the dominant arithmetic cost.  Expanding the component products gives
%   sum_{s,t} [(conj(L_s).*w_v)*L_t.'] .* [(R_s.*w_c)*R_t'].
% This uses only Nmu-by-Nmu intermediates and BLAS GEMMs.
nmu = size(left{1}, 1);
value = complex(zeros(nmu, nmu, 'like', left{1}));
occ_weight = reshape(occ_weight, 1, []);
unocc_weight = reshape(unocc_weight, 1, []);
for ileft_component = 1:numel(left)
    weighted_left = bsxfun(@times, conj(left{ileft_component}), ...
        occ_weight);
    weighted_right = bsxfun(@times, right{ileft_component}, ...
        unocc_weight);
    for iright_component = 1:numel(left)
        left_gram = weighted_left * left{iright_component}.';
        right_gram = weighted_right * right{iright_component}';
        value = value + left_gram .* right_gram;
    end
end
end

function [center, radius, ok] = local_contour(ev_occ, ev_unocc)
center = 0.5 * (ev_occ(1) + ev_occ(end));
half_width = 0.5 * (ev_occ(end) - ev_occ(1));
max_radius = min(abs(ev_unocc - center));
ok = isfinite(max_radius) && max_radius > half_width;
if ~ok
    radius = NaN;
    return;
end
radius = 0.5 * (half_width + max_radius);
if radius >= max_radius
    radius = 0.9 * max_radius;
end
end
