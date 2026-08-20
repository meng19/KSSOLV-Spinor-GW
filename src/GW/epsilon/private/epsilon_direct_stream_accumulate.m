function acc = epsilon_direct_stream_accumulate(ctx, acc, block)
%EPSILON_DIRECT_STREAM_ACCUMULATE Direct save_mem accumulation.

for iv_local = 1:numel(block.valence_bands)
    iv = block.valence_bands(iv_local);
    for ic_local = 1:numel(block.conduction_bands)
        ic = block.conduction_bands(ic_local);
        vector = getm_epsilon(iv, ic, block.wfnkq, block.wfnk, ...
            block.fft, block.idx, block.ispin, ctx.nspinor, ctx.use_gpu);
        eden_pages = zeros(1, 1, ctx.pol.nfreq);
        for ifreq = 1:ctx.pol.nfreq
            freq = ctx.pol.freq(ifreq) / ctx.ryd;
            eden_pages(1, 1, ifreq) = get_eden( ...
                iv, ic, block.wfnkq, block.wfnk, block.ispin, ...
                ctx.options, freq, ctx.eps);
        end
        acc.chi0 = epsilon_add_state_batch( ...
            acc.chi0, vector, eden_pages, block.g_maps);
    end
end
end
