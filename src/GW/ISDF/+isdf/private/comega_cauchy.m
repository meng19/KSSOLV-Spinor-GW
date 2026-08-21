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
products = component_products(left, right, [], []);
previous = [];
relative_error = inf;
for iteration = 1:options.MaxIter
    npoints = 2^(iteration + 3);
    result = complex(zeros(nmu, nmu, 'like', left{1}));
    for ipoint = 0:npoints-1
        theta = 2 * pi * ipoint / npoints;
        exp_theta = exp(1i * theta);
        z = center + radius * exp_theta;
        occ_weight = 1 ./ (z - ev_occ_work);
        unocc_weight = 1 ./ (z - ev_unocc_work);
        result = result + local_weighted_products( ...
            products, occ_weight, unocc_weight) * ...
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

function value = local_weighted_products(products, occ_weight, unocc_weight)
pair_weight = reshape(occ_weight(:) * unocc_weight(:).', 1, []);
weighted_products = bsxfun(@times, products, pair_weight);
value = weighted_products * products';
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
