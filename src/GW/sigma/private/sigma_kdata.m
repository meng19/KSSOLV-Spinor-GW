function kdata = sigma_kdata(ctx, ik)
%SIGMA_KDATA Build irreducible q data for one sigma k-point.

kdata.rk = ctx.sig.qpt(ik, :);
kdata.syms = subgrp(kdata.rk, ctx.syms);
[kdata.nrk, kdata.neq, kdata.indrk] = irrbz(kdata.syms, ctx.gr);
if ctx.sig.no_symmetries_q_grid
    kdata.nrk = ctx.gr.nf;
    kdata.indrk = 1:kdata.nrk;
    kdata.neq = ones(1, kdata.nrk);
end
end
