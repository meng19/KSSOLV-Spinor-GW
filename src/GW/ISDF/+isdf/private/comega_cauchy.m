function [result, relative_error, iteration, fallback_direct] = comega_cauchy( ...
    left, right, ev_occ, ev_unocc, options)
%COMEGA_CAUCHY Cauchy reduced polarizability for frequency pages.
%
% The separable node weight w(v,c; z) = wv(v; z)*wc(c; z) lets every
% quadrature node of a Newton iteration be folded into one effective pair
% weight gamma(v,c) = sum_z c_z*wv(v;z)*wc(c;z).  Each iteration then costs
% one pair of GEMMs on prebuilt coefficient matrices instead of two GEMMs
% per node:
%   result = sum_b (P_b .* gamma_b) * Q_b'
%   P_b(p; v,c) = sum_s conj(L_s(p,v)).*R_s(p,c)
%   Q_b(q; v,c) = sum_t L_t(q,v).*R_t(q,c)
% P and Q depend only on the block coefficients, not on frequency or node,
% so they are built once per call and reused by the static, plus-, and
% minus-frequency pages.  Column chunks over the conduction index bound the
% working-set memory on large product spaces.

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
chunks = local_coefficient_chunks(left, right);
for ifreq = 1:nfreq
    omega = freq(ifreq);
    if omega == 0
        [result(:, :, ifreq), relative_errors(ifreq), iterations(ifreq)] = ...
            local_cauchy_resolvent(chunks, ev_occ, ev_unocc, options);
        continue;
    end

    [plus_page, plus_error, plus_iter, plus_ok] = ...
        local_cauchy_resolvent( ...
        chunks, ev_occ, ev_unocc - omega, options);
    [minus_page, minus_error, minus_iter, minus_ok] = ...
        local_cauchy_resolvent( ...
        chunks, ev_occ, ev_unocc + omega, options);
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

function chunks = local_coefficient_chunks(left, right)
% Fold the component sum into nmu-by-(Nv*Nc) coefficient matrices, chunked
% over the conduction index so the stored working set stays bounded.
nmu = size(left{1}, 1);
nv = size(left{1}, 2);
nc = size(right{1}, 2);
% Each stored column pair costs 2 * 16 * nmu bytes (one P and one Q column).
max_chunk_columns = max(1, floor(5.12e8 / (32 * nmu)));
jchunk = max(1, floor(max_chunk_columns / nv));
chunks = struct('P', {}, 'Q', {}, 'first', {}, 'last', {});
for jfirst = 1:jchunk:nc
    jlast = min(nc, jfirst + jchunk - 1);
    P = complex(zeros(nmu, nv * (jlast - jfirst + 1), 'like', left{1}));
    Q = complex(zeros(nmu, nv * (jlast - jfirst + 1), 'like', right{1}));
    for ic = jfirst:jlast
        local_columns = (ic - jfirst) * nv + (1:nv);
        P_slice = complex(zeros(nmu, nv, 'like', left{1}));
        Q_slice = complex(zeros(nmu, nv, 'like', right{1}));
        for icomponent = 1:numel(left)
            P_slice = P_slice + bsxfun(@times, conj(left{icomponent}), ...
                right{icomponent}(:, ic));
            Q_slice = Q_slice + bsxfun(@times, left{icomponent}, ...
                conj(right{icomponent}(:, ic)));
        end
        P(:, local_columns) = P_slice;
        Q(:, local_columns) = Q_slice;
    end
    chunks(end + 1) = struct('P', {P}, 'Q', {Q}, ...
        'first', (jfirst - 1) * nv + 1, 'last', jlast * nv); %#ok<AGROW>
end
end

function [result, relative_error, iteration, ok] = local_cauchy_resolvent( ...
    chunks, ev_occ, ev_unocc, options)
[center, radius, ok] = local_contour(ev_occ, ev_unocc);
nmu = size(chunks(1).P, 1);
if ~ok
    result = complex(zeros(nmu, nmu, 'like', chunks(1).P));
    relative_error = inf;
    iteration = 0;
    return;
end
ev_occ_work = ev_occ;
ev_unocc_work = ev_unocc;
if isa(chunks(1).P, 'gpuArray')
    ev_occ_work = gpuArray(ev_occ_work);
    ev_unocc_work = gpuArray(ev_unocc_work);
end
previous = [];
gamma_previous = [];
relative_error = inf;
for iteration = 1:options.MaxIter
    npoints = 2^(iteration + 3);
    if isempty(previous)
        point_indices = 0:npoints-1;
        gamma = zeros(numel(ev_occ), numel(ev_unocc));
    else
        % The nodes for npoints/2 are the even nodes for npoints.  Their
        % completed trapezoidal contribution is folded in through the pair
        % weights, so only the new odd nodes are evaluated.
        point_indices = 1:2:npoints-1;
        gamma = 0.5 * gamma_previous;
    end
    % Node separability stacks the new node weights into one small
    % Nv-by-Nc accumulation instead of a per-node GEMM sweep.
    exp_theta = exp(1i * (2 * pi * point_indices / npoints));
    z = center + radius * exp_theta;
    occ_weight = 1 ./ (z(:) - ev_occ_work(:).');       % nnodes-by-Nv
    unocc_weight = 1 ./ (z(:) - ev_unocc_work(:).');   % nnodes-by-Nc
    gamma = gamma + occ_weight.' * ...
        (unocc_weight .* (radius * exp_theta(:) / npoints));
    result = local_folded_apply(chunks, gamma);
    if ~isempty(previous)
        relative_error = gather_if_gpu( ...
            norm(result - previous, 'fro') / max(1, norm(result, 'fro')));
        if relative_error <= options.froErr
            return;
        end
    end
    previous = result;
    gamma_previous = gamma;
end
end

function value = local_folded_apply(chunks, gamma)
% result = P*diag(gamma(:))*Q' accumulated over the stored chunks.  The
% column ranges of each chunk match the column-major ordering of gamma.
gamma_row = reshape(gamma, 1, []);
value = complex(zeros(size(chunks(1).P, 1), size(chunks(1).P, 1), ...
    'like', chunks(1).P));
for ichunk = 1:numel(chunks)
    scaled = chunks(ichunk).P .* gamma_row( ...
        chunks(ichunk).first:chunks(ichunk).last);
    value = value + scaled * chunks(ichunk).Q.';
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
