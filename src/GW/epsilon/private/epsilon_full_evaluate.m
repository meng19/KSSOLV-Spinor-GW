function contribution = epsilon_full_evaluate(ctx, block, use_isdf)
%EPSILON_FULL_EVALUATE Build full matrix-element contributions.

if ctx.save_mem && ~use_isdf
    contribution.stream_direct = true;
    return;
end

gme_size = [numel(block.idx.q), numel(block.valence_bands), ...
    numel(block.conduction_bands)];
npairs = max(1, numel(block.valence_bands) * ...
    numel(block.conduction_bands));
progress_work = local_progress_work(block);
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
        epsilon_progress(block, ...
            progress_work * 0.25 * ispinor / ctx.nspinor, ...
            sprintf('Eps q:%d ik:%d s:%d', ...
            block.iq, block.ik, block.ispin));
    end
    isdf_options = ctx.eps.isdf;
    if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
        isdf_options.rank = ceil(sqrt(npairs) * ctx.eps.isdf.rank_ratio);
    end
    contribution.gme = isdf.matrix_elements(left, right, ...
        block.idx.q, size(block.fft.Nfft1), isdf_options);
    epsilon_progress(block, progress_work * 0.65, ...
        sprintf('Eps q:%d ik:%d s:%d', block.iq, block.ik, block.ispin));
else
    pair_count = 0;
    for iv_local = 1:numel(block.valence_bands)
        iv = block.valence_bands(iv_local);
        for ic_local = 1:numel(block.conduction_bands)
            ic = block.conduction_bands(ic_local);
            contribution.gme(:, iv_local, ic_local) = getm_epsilon( ...
                iv, ic, block.wfnkq, block.wfnk, block.fft, block.idx, ...
                block.ispin, ctx.nspinor, ctx.use_gpu);
            pair_count = pair_count + 1;
            epsilon_progress(block, ...
                progress_work * 0.5 * pair_count / npairs, ...
                sprintf('Eps q:%d ik:%d s:%d v:%d c:%d', ...
                block.iq, block.ik, block.ispin, iv, ic));
        end
    end
end

contribution.eden = zeros(numel(block.valence_bands), ...
    numel(block.conduction_bands), ctx.pol.nfreq);
pair_count = 0;
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
        pair_count = pair_count + 1;
        epsilon_progress(block, ...
            progress_work * (0.5 + 0.5 * pair_count / npairs), ...
            sprintf('Eps q:%d ik:%d s:%d v:%d c:%d', ...
            block.iq, block.ik, block.ispin, iv, ic));
    end
end
end

function work = local_progress_work(block)
if isfield(block, 'progress') && isfield(block.progress, 'block_work')
    work = block.progress.block_work;
else
    work = max(1, numel(block.valence_bands) * ...
        numel(block.conduction_bands));
end
end
