function block = epsilon_prepare_block(ctx, iq, ik, ispin, prepared)
qdata = ctx.qdata{iq};
qq = ctx.pol.qpt(iq, :);
rk = ctx.gr.f(qdata.indrk(ik), :);

if nargin >= 5 && ~isempty(prepared)
    wfnk = prepared.wfnk;
    wfnkq = prepared.wfnkq;
    fft = prepared.fft;
    idx = prepared.idx;
else
    wfnk = genwf(rk, ctx.gr, ctx.gvec, ctx.syms, ctx.sys, ...
        ctx.options, ctx.wfc_cutoff, ctx.nbands, ctx.use_gpu);
    wfnkq = genwf(rk + qq, ctx.gr, ctx.gvec, ctx.syms, ctx.sys, ...
        ctx.options, ctx.wfc_cutoff, ctx.nbands, ctx.use_gpu);
    [fft, idx] = epsilon_prefft( ...
        wfnkq, wfnk, iq, ik, ctx.pol, [], [], ctx.use_gpu);
end

block.iq = iq;
block.ik = ik;
block.ispin = ispin;
block.ik_fbz = qdata.indrk(ik);
block.q = qq;
block.k = rk;
block.weight = qdata.neq(ik);
block.g_maps = qdata.g_maps{ik};
block.star_kpoints = qdata.rqs{ik};
block.wfnk = wfnk;
block.wfnkq = wfnkq;
block.fft = fft;
block.idx = idx;
block.occ_vkq = get_occ(ctx.options, wfnkq.ikq, ispin);
block.occ_ck = get_occ(ctx.options, wfnk.ikq, ispin);
block.valence_bands = 1:sum(block.occ_vkq > 0);
block.conduction_bands = (sum(block.occ_ck > 0) + 1):ctx.nbands;
block.ev_occ = ctx.options.ev( ...
    block.valence_bands, wfnkq.ikq, ispin);
block.ev_unocc = ctx.options.ev( ...
    block.conduction_bands, wfnk.ikq, ispin);
end
