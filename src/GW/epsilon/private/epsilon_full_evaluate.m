function contribution = epsilon_full_evaluate(ctx, block, use_isdf)
%EPSILON_FULL_EVALUATE Build full matrix-element contributions.

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
