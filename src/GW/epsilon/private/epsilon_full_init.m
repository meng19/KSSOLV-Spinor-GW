function acc = epsilon_full_init(ctx, iq)
%EPSILON_FULL_INIT Initialize full G-space epsilon accumulation.

nmtx = ctx.pol.nmtx(iq);
if ctx.use_gpu
    acc.chi0 = gpuArray.zeros(nmtx, nmtx, ctx.pol.nfreq);
else
    acc.chi0 = zeros(nmtx, nmtx, ctx.pol.nfreq);
end
acc.deferred = cell(ctx.nspin, ctx.gr.nf);
end
