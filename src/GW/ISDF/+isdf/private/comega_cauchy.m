function [result, relative_error, iteration] = comega_cauchy( ...
    left, right, ev_occ, ev_unocc, options)
%COMEGA_CAUCHY Static Cauchy reduced polarizability.

gap_min = ev_unocc(1) - ev_occ(end);
if gap_min <= 0
    error('ISDF:CauchyNoGap', ...
        'Cauchy integral requires a positive band gap.');
end
center = 0.5 * (ev_occ(1) + ev_occ(end));
half_width = 0.5 * (ev_occ(end) - ev_occ(1));
radius = half_width + 0.5 * gap_min;
max_radius = ev_unocc(1) - center;
if radius >= max_radius
    radius = 0.9 * max_radius;
end

nmu = size(left{1}, 1);
ncomponents = numel(left);
if ncomponents ~= numel(right)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
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
    result = complex(zeros(nmu, nmu, 'like', left{1}));
    for ipoint = 0:npoints-1
        theta = 2 * pi * ipoint / npoints;
        exp_theta = exp(1i * theta);
        z = center + radius * exp_theta;
        for icomponent = 1:ncomponents
            left_s = left{icomponent};
            right_s = right{icomponent};
            for jcomponent = 1:ncomponents
                left_t = left{jcomponent};
                right_t = right{jcomponent};
                occ_matrix = conj(left_s) * ...
                    diag(1 ./ (z - ev_occ_work)) * left_t.';
                unocc_matrix = right_s * ...
                    diag(1 ./ (z - ev_unocc_work)) * right_t';
                result = result + (occ_matrix .* unocc_matrix) * ...
                    (radius * exp_theta / npoints);
            end
        end
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
