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
        error('ISDF:ReducedEpsilonNotIntegrated', ...
            'Reduced epsilon is integrated in the next migration task.');
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

function eps = local_finalize_full(ctx, eps, acc, iq)
chi0 = acc.chi0;
if ~ctx.save_mem
    for ispin = 1:size(acc.deferred, 1)
        for ik_fbz = 1:size(acc.deferred, 2)
            entries = acc.deferred{ispin, ik_fbz};
            for ientry = 1:numel(entries)
                entry = entries{ientry};
                chi0 = local_add_states(chi0, entry{1}, entry{2});
            end
        end
    end
end

chi0 = chi0 * ctx.fact;
nmtx = ctx.pol.nmtx(iq);
coulg = coulG_select(eps, nmtx, ctx.pol.isrtx(:, iq), ...
    ctx.ekin(:, iq), 0, ctx.pol.mtx{:, iq}, ctx.gvec, ctx.sys, iq);

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
