function acc = epsilon_reduced_init( ...
    ctx, iq, need_full_inverse, need_screened_w)
%EPSILON_REDUCED_INIT Initialize reduced-basis epsilon accumulation.

acc.need_full_inverse = need_full_inverse;
acc.need_screened_w = need_screened_w;
if acc.need_full_inverse
    if ctx.use_gpu
        acc.chi0 = gpuArray.zeros(ctx.pol.nmtx(iq), ...
            ctx.pol.nmtx(iq), ctx.pol.nfreq);
    else
        acc.chi0 = zeros(ctx.pol.nmtx(iq), ctx.pol.nmtx(iq), ...
            ctx.pol.nfreq);
    end
else
    acc.chi0 = [];
end
acc.zeta_blocks = {};
acc.coeff_blocks = {};
acc.rank = cell(ctx.nspin, ctx.qdata{iq}.nrq);
acc.info = cell(ctx.nspin, ctx.qdata{iq}.nrq);
end
