function ops = epsilon_ops(ctx)
ops.name = ctx.method;
switch ctx.method
    case 'direct'
        ops.init = @(iq) local_init_full(ctx, iq);
        ops.evaluate = @(block) local_evaluate_full(ctx, block, false);
        ops.accumulate = @(acc, contribution, block) ...
            local_accumulate_full(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) local_finalize_full(ctx, eps, acc, iq);
    case 'matrix_elements'
        ops.init = @(iq) local_init_full(ctx, iq);
        ops.evaluate = @(block) local_evaluate_full(ctx, block, true);
        ops.accumulate = @(acc, contribution, block) ...
            local_accumulate_full(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) local_finalize_full(ctx, eps, acc, iq);
    case 'reduced_basis'
        output_mode = lower(ctx.eps.isdf.output);
        need_full_inverse = any(strcmp( ...
            output_mode, {'full_inverse', 'both'}));
        need_screened_w = any(strcmp( ...
            output_mode, {'screened_w', 'both'}));
        if ~need_full_inverse && ~need_screened_w
            error('ISDF:ReducedEpsilonOutput', ...
                'Unknown ISDF reduced-basis epsilon output "%s".', ...
                ctx.eps.isdf.output);
        end
        ops.init = @(iq) local_init_reduced( ...
            ctx, iq, need_full_inverse, need_screened_w);
        ops.evaluate = @(block) local_evaluate_reduced(ctx, block);
        ops.accumulate = @(acc, contribution, block) ...
            local_accumulate_reduced(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) ...
            local_finalize_reduced(ctx, eps, acc, iq);
end
end

function acc = local_init_reduced( ...
    ctx, iq, need_full_inverse, need_screened_w)
acc.need_full_inverse = need_full_inverse;
acc.need_screened_w = need_screened_w;
if acc.need_full_inverse
    acc.chi0 = zeros(ctx.pol.nmtx(iq));
else
    acc.chi0 = [];
end
acc.zeta_blocks = {};
acc.coeff_blocks = {};
acc.rank = cell(ctx.nspin, ctx.qdata{iq}.nrq);
acc.info = cell(ctx.nspin, ctx.qdata{iq}.nrq);
end

function contribution = local_evaluate_reduced(ctx, block)
left = cell(1, ctx.nspinor);
right = cell(1, ctx.nspinor);
for ispinor = 1:ctx.nspinor
    left{ispinor} = isdf.real_component(block.wfnkq, ...
        block.fft.Nfft1, block.idx.kq, block.ispin, ispinor, ...
        block.valence_bands);
    right{ispinor} = isdf.real_component(block.wfnk, ...
        block.fft.Nfft2, block.idx.k, block.ispin, ispinor, ...
        block.conduction_bands);
end
isdf_options = ctx.eps.isdf;
if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
    npairs = numel(block.valence_bands) * ...
        numel(block.conduction_bands);
    isdf_options.rank = ceil(sqrt(npairs) * ctx.eps.isdf.rank_ratio);
end
space = isdf.build_space(left, right, block.idx.q, ...
    size(block.fft.Nfft1), isdf_options);
solver.method = ctx.eps.isdf.reduced_solver;
solver.froErr = ctx.eps.isdf.cauchy_froErr;
solver.MaxIter = ctx.eps.isdf.cauchy_MaxIter;
polar = isdf.polarizability( ...
    space, block.ev_occ, block.ev_unocc, solver);
contribution.space = space;
contribution.polar = polar;
end

function acc = local_accumulate_reduced(~, acc, contribution, block)
for it = 1:numel(block.g_maps)
    zeta_star = contribution.space.zeta_g(block.g_maps{it}, :);
    zeta_chi = conj(zeta_star);
    coeff_chi = conj(contribution.polar.coeff);
    if acc.need_full_inverse
        acc.chi0 = acc.chi0 + zeta_chi * coeff_chi * zeta_chi';
    end
    if acc.need_screened_w
        acc.zeta_blocks{end + 1} = zeta_chi;
        acc.coeff_blocks{end + 1} = coeff_chi;
    end
end
acc.rank{block.ispin, block.ik} = contribution.space.rank;
acc.info{block.ispin, block.ik} = contribution.polar.info;
end

function eps = local_finalize_reduced(ctx, eps, acc, iq)
for ispin = 1:size(acc.rank, 1)
    for ik = 1:size(acc.rank, 2)
        eps.isdf_reduced_rank{iq, ispin, ik} = acc.rank{ispin, ik};
        eps.isdf_reduced_info{iq, ispin, ik} = acc.info{ispin, ik};
    end
end

coulg = coulG_select(ctx.eps, ctx.pol.nmtx(iq), ...
    ctx.pol.isrtx(:, iq), ctx.ekin(:, iq), 0, ...
    ctx.pol.mtx{:, iq}, ctx.gvec, ctx.sys, iq);
if acc.need_full_inverse
    if iq == 1
        eps.inv = cell(ctx.nq, 1);
    end
    epsilon_matrix = eye(ctx.pol.nmtx(iq)) - ...
        coulg(:) .* (ctx.fact * acc.chi0);
    eps.inv{iq} = inv(epsilon_matrix);
end
if acc.need_screened_w && ~isfield(eps, 'isdf_screened_w')
    eps.isdf_screened_w = cell(ctx.nq, 1);
end
if acc.need_screened_w && ~isempty(acc.zeta_blocks)
    combined_space.zeta_g = cat(2, acc.zeta_blocks{:});
    combined_polar.coeff = ctx.fact * blkdiag(acc.coeff_blocks{:});
    eps.isdf_screened_w{iq} = isdf.screened_w( ...
        combined_space, coulg(:), combined_polar);
end
end

function acc = local_init_full(ctx, iq)
nmtx = ctx.pol.nmtx(iq);
if ctx.use_gpu
    acc.chi0 = gpuArray.zeros(nmtx, nmtx, ctx.pol.nfreq);
else
    acc.chi0 = zeros(nmtx, nmtx, ctx.pol.nfreq);
end
acc.deferred = cell(ctx.nspin, ctx.gr.nf);
end

function contribution = local_evaluate_full(ctx, block, use_isdf)
if ctx.save_mem && ~use_isdf
    contribution.stream_direct = true;
    return;
end

gme_size = [numel(block.idx.q), numel(block.valence_bands), ...
    numel(block.conduction_bands)];
if ctx.use_gpu
    contribution.gme = gpuArray.zeros(gme_size);
else
    contribution.gme = zeros(gme_size);
end

if use_isdf
    left = cell(1, ctx.nspinor);
    right = cell(1, ctx.nspinor);
    for ispinor = 1:ctx.nspinor
        left{ispinor} = isdf.real_component(block.wfnkq, ...
            block.fft.Nfft1, block.idx.kq, block.ispin, ispinor, ...
            block.valence_bands);
        right{ispinor} = isdf.real_component(block.wfnk, ...
            block.fft.Nfft2, block.idx.k, block.ispin, ispinor, ...
            block.conduction_bands);
    end
    isdf_options = ctx.eps.isdf;
    if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
        isdf_options.rank = ceil(sqrt(numel(block.valence_bands) * ...
            numel(block.conduction_bands)) * ctx.eps.isdf.rank_ratio);
    end
    contribution.gme = isdf.matrix_elements(left, right, ...
        block.idx.q, size(block.fft.Nfft1), isdf_options);
else
    for iv_local = 1:numel(block.valence_bands)
        iv = block.valence_bands(iv_local);
        for ic_local = 1:numel(block.conduction_bands)
            ic = block.conduction_bands(ic_local);
            contribution.gme(:, iv_local, ic_local) = getm_epsilon( ...
                iv, ic, block.wfnkq, block.wfnk, block.fft, block.idx, ...
                block.ispin, ctx.nspinor, ctx.use_gpu);
        end
    end
end

contribution.eden = zeros(numel(block.valence_bands), ...
    numel(block.conduction_bands), ctx.pol.nfreq);
for iv_local = 1:numel(block.valence_bands)
    iv = block.valence_bands(iv_local);
    for ic_local = 1:numel(block.conduction_bands)
        ic = block.conduction_bands(ic_local);
        for ifreq = 1:ctx.pol.nfreq
            freq = ctx.pol.freq(ifreq) / ctx.ryd;
            contribution.eden(iv_local, ic_local, ifreq) = get_eden( ...
                iv, ic, block.wfnkq, block.wfnk, block.ispin, ...
                ctx.options, freq, ctx.eps);
        end
    end
end
end

function acc = local_accumulate_full(ctx, acc, contribution, block)
if isfield(contribution, 'stream_direct')
    acc = local_accumulate_direct_stream(ctx, acc, block);
    return;
end

for it = 1:numel(block.g_maps)
    gme = contribution.gme(block.g_maps{it}, :, :);
    if ctx.save_mem
        acc.chi0 = local_add_states(acc.chi0, gme, contribution.eden);
    elseif isempty(acc.deferred{block.ispin, block.ik_fbz})
        acc.deferred{block.ispin, block.ik_fbz} = ...
            {{gme, contribution.eden}};
    else
        acc.deferred{block.ispin, block.ik_fbz}{end + 1} = ...
            {gme, contribution.eden};
    end
end
end

function acc = local_accumulate_direct_stream(ctx, acc, block)
for iv_local = 1:numel(block.valence_bands)
    iv = block.valence_bands(iv_local);
    for ic_local = 1:numel(block.conduction_bands)
        ic = block.conduction_bands(ic_local);
        vector = getm_epsilon(iv, ic, block.wfnkq, block.wfnk, ...
            block.fft, block.idx, block.ispin, ctx.nspinor, ctx.use_gpu);
        for ifreq = 1:ctx.pol.nfreq
            freq = ctx.pol.freq(ifreq) / ctx.ryd;
            eden = get_eden(iv, ic, block.wfnkq, block.wfnk, ...
                block.ispin, ctx.options, freq, ctx.eps);
            for it = 1:numel(block.g_maps)
                gme = vector(block.g_maps{it});
                acc.chi0(:, :, ifreq) = acc.chi0(:, :, ifreq) + ...
                    conj(gme) * gme.' * eden;
            end
        end
    end
end
end

function eps = local_finalize_full(ctx, eps, acc, iq)
chi0 = acc.chi0;
if ~ctx.save_mem
    for ispin = 1:size(acc.deferred, 1)
        for ik_fbz = 1:size(acc.deferred, 2)
            entries = acc.deferred{ispin, ik_fbz};
            for ientry = 1:numel(entries)
                entry = entries{ientry};
                chi0 = local_add_deferred_states(chi0, entry{1}, entry{2});
            end
        end
    end
end

chi0 = chi0 * ctx.fact;
nmtx = ctx.pol.nmtx(iq);
coulg = coulG_select(eps, nmtx, ctx.pol.isrtx(:, iq), ...
    ctx.ekin(:, iq), 0, ctx.pol.mtx{:, iq}, ctx.gvec, ctx.sys, iq);

if iq == 1
    eps.inv = cell(ctx.nq, 1);
end

if ctx.use_gpu
    if ~isa(coulg, 'gpuArray')
        coulg_gpu = gpuArray(coulg(:));
    else
        coulg_gpu = coulg(:);
    end
    identity = eye(nmtx, 'gpuArray');
    epsilon_pages = repmat(identity, 1, 1, ctx.pol.nfreq);
    epsilon_pages = epsilon_pages - bsxfun(@times, coulg_gpu, chi0);
    eps.inv{iq} = gather(pagefun(@inv, epsilon_pages));
else
    identity = repmat(eye(nmtx), 1, 1, ctx.pol.nfreq);
    epsilon_pages = identity - bsxfun(@times, coulg(:), chi0);
    inverse_pages = cell(1, size(epsilon_pages, 3));
    for ifreq = 1:size(epsilon_pages, 3)
        inverse_pages{ifreq} = inv(epsilon_pages(:, :, ifreq));
    end
    eps.inv{iq} = cat(3, inverse_pages{:});
end
end

function chi0 = local_add_states(chi0, gme, eden)
for iv = 1:size(gme, 2)
    for ic = 1:size(gme, 3)
        vector = gme(:, iv, ic);
        for ifreq = 1:size(eden, 3)
            chi0(:, :, ifreq) = chi0(:, :, ifreq) + ...
                conj(vector) * vector.' * eden(iv, ic, ifreq);
        end
    end
end
end

function chi0 = local_add_deferred_states(chi0, gme, eden)
gme = reshape(gme, size(gme, 1), []);
eden = reshape(eden, [], size(eden, 3));
for ifreq = 1:size(eden, 2)
    weight = eden(:, ifreq).';
    if isa(gme, 'gpuArray') && ~isa(weight, 'gpuArray')
        weight = gpuArray(weight);
    end
    weighted_gme = bsxfun(@times, gme, weight);
    chi0(:, :, ifreq) = chi0(:, :, ifreq) + ...
        conj(gme) * weighted_gme.';
end
end
