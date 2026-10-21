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
freq = options.freq(:).';
nmu = size(left{1}, 1);
nfreq = numel(freq);
result = complex(zeros(nmu, nmu, nfreq, 'like', left{1}));
relative_errors = zeros(1, nfreq);
iterations = zeros(1, nfreq);
fallback_direct = false(1, nfreq);
chunks = local_coefficient_chunks(left, right, options);
for ifreq = 1:nfreq
    omega = freq(ifreq);
    if omega == 0
        [page, page_error, page_iter, page_ok] = ...
            local_cauchy_resolvent(chunks, left, right, ev_occ, ev_unocc, options);
        if page_ok
            result(:, :, ifreq) = page;
            relative_errors(ifreq) = page_error;
            iterations(ifreq) = page_iter;
        else
            % A zero/negative gap cannot be enclosed by this contour.  Do
            % not silently return the zero placeholder from the failed
            % contour construction; retain the exact reduced summation.
            result(:, :, ifreq) = local_direct_resolvent( ...
                chunks, left, right, ev_occ, ev_unocc);
            fallback_direct(ifreq) = true;
        end
        continue;
    end

    [plus_page, plus_error, plus_iter, plus_ok] = ...
        local_cauchy_resolvent( ...
        chunks, left, right, ev_occ, ev_unocc - omega, options);
    [minus_page, minus_error, minus_iter, minus_ok] = ...
        local_cauchy_resolvent( ...
        chunks, left, right, ev_occ, ev_unocc + omega, options);
    if plus_ok && minus_ok
        result(:, :, ifreq) = 0.5 * (plus_page + minus_page);
        relative_errors(ifreq) = max(plus_error, minus_error);
        iterations(ifreq) = max(plus_iter, minus_iter);
    else
        % Keep a successfully converged half of the resolvent identity.
        % For real positive frequencies the minus contour commonly remains
        % valid while the plus contour crosses a resonance, so replacing
        % only the failed half avoids an unnecessary exact contraction.
        if plus_ok
            plus_page_exact = plus_page;
            page_plus_error = plus_error;
            page_plus_iter = plus_iter;
        else
            plus_page_exact = local_direct_resolvent( ...
                chunks, left, right, ev_occ, ev_unocc - omega);
            page_plus_error = 0;
            page_plus_iter = 0;
        end
        if minus_ok
            minus_page_exact = minus_page;
            page_minus_error = minus_error;
            page_minus_iter = minus_iter;
        else
            minus_page_exact = local_direct_resolvent( ...
                chunks, left, right, ev_occ, ev_unocc + omega);
            page_minus_error = 0;
            page_minus_iter = 0;
        end
        result(:, :, ifreq) = 0.5 * (plus_page_exact + minus_page_exact);
        relative_errors(ifreq) = max(page_plus_error, page_minus_error);
        iterations(ifreq) = max(page_plus_iter, page_minus_iter);
        fallback_direct(ifreq) = true;
    end
end
relative_error = max(relative_errors);
iteration = max(iterations);
end

function chunks = local_coefficient_chunks(left, right, options)
% Describe pair-column chunks; cache their P/Q values only when small.
nmu = size(left{1}, 1);
npairs = size(left{1}, 2) * size(right{1}, 2);
% One pair column stores P and Q, i.e. 2 complex-double vectors of nmu.
pair_block = max(1, floor(options.cauchy_coefficient_block_bytes / ...
    max(1, 32 * nmu)));
pair_block = min(npairs, pair_block);
cache_all = local_cache_coefficients(options, nmu, npairs);
chunks = struct('P', {}, 'Q', {}, 'first', {}, 'last', {});
for first = 1:pair_block:npairs
    last = min(npairs, first + pair_block - 1);
    chunk = struct('P', [], 'Q', [], 'first', first, 'last', last);
    if cache_all
        [chunk.P, chunk.Q] = local_chunk_coefficients(chunk, left, right);
    end
    chunks(end + 1) = chunk; %#ok<AGROW>
end
end

function tf = local_cache_coefficients(options, nmu, npairs)
switch options.cauchy_coefficient_storage
    case 'full'
        tf = true;
    case 'stream'
        tf = false;
    otherwise
        tf = 32 * nmu * npairs <= options.cauchy_coefficient_cache_bytes;
end
end

function [result, relative_error, iteration, ok] = local_cauchy_resolvent( ...
    chunks, left, right, ev_occ, ev_unocc, options)
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
gamma_previous = [];
relative_error = inf;
% ``ok`` means both a valid contour and successful quadrature convergence.
% Leaving it true after MaxIter would silently accept an inaccurate page.
ok = false;
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
    result = local_folded_apply(chunks, gamma, left, right);
    if ~isempty(previous)
        relative_error = gather_if_gpu( ...
            norm(result - previous, 'fro') / max(1, norm(result, 'fro')));
        if relative_error <= options.froErr
            ok = true;
            return;
        end
    end
    previous = result;
    gamma_previous = gamma;
end
end

function value = local_direct_resolvent(chunks, left, right, ev_occ, ev_unocc)
% Exact reduced pair summation, evaluated through the cached P/Q chunks.
% This is algebraically identical to comega_direct for one resolvent page,
% but avoids reconstructing one real-space product at a time.
ev_occ_work = ev_occ;
ev_unocc_work = ev_unocc;
if isa(left{1}, 'gpuArray')
    ev_occ_work = gpuArray(ev_occ_work);
    ev_unocc_work = gpuArray(ev_unocc_work);
end
gamma = 1 ./ (ev_occ_work(:) - ev_unocc_work(:).');
value = local_folded_apply(chunks, gamma, left, right);
end

function value = local_folded_apply(chunks, gamma, left, right)
% Accumulate P*diag(gamma(:))*Q'; rebuild streamed chunks one at a time.
gamma_row = reshape(gamma, 1, []);
value = complex(zeros(size(left{1}, 1), size(left{1}, 1), ...
    'like', left{1}));
for ichunk = 1:numel(chunks)
    [P, Q] = local_chunk_coefficients(chunks(ichunk), left, right);
    scaled = P .* gamma_row( ...
        chunks(ichunk).first:chunks(ichunk).last);
    value = value + scaled * Q.';
end
end

function [P, Q] = local_chunk_coefficients(chunk, left, right)
if ~isempty(chunk.P)
    P = chunk.P;
    Q = chunk.Q;
    return;
end
nmu = size(left{1}, 1);
nv = size(left{1}, 2);
columns = chunk.first:chunk.last;
iv = mod(columns - 1, nv) + 1;
ic = floor((columns - 1) / nv) + 1;
P = complex(zeros(nmu, numel(columns), 'like', left{1}));
Q = complex(zeros(nmu, numel(columns), 'like', right{1}));
for icomponent = 1:numel(left)
    P = P + conj(left{icomponent}(:, iv)) .* right{icomponent}(:, ic);
    Q = Q + left{icomponent}(:, iv) .* conj(right{icomponent}(:, ic));
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
