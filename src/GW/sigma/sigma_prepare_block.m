function block = sigma_prepare_block(ctx, ik, iq, in, ispin, prepared)
%SIGMA_PREPARE_BLOCK Normalize one target-band/full-BZ-q sigma block.

kdata = ctx.kdata{ik};
iq_fbz = kdata.indrk(iq);
rk = kdata.rk;
qq = ctx.gr.f(iq_fbz, :);

if nargin < 6 || isempty(prepared)
    prepared = struct();
end
if isfield(prepared, 'wfnk')
    wfnk = prepared.wfnk;
else
    wfnk = genwf(rk, ctx.gr, ctx.gvec, ctx.syms, ctx.sys, ...
        ctx.options, ctx.wfc_cutoff, ctx.nbands, ctx.use_gpu);
end
if isfield(prepared, 'wfnkq')
    wfnkq = prepared.wfnkq;
else
    wfnkq = genwf(rk - qq, ctx.gr, ctx.gvec, ctx.syms, ctx.sys, ...
        ctx.options, ctx.wfc_cutoff, ctx.nbands, ctx.use_gpu);
end
if isfield(prepared, 'idx')
    idx = prepared.idx;
else
    idx = sigma_prefft(wfnkq, wfnk, ctx.fbz.mtx{:, iq_fbz}, ...
        iq, ik, ctx.sys, [], ctx.use_gpu);
end

coulg = coulG_select(ctx.sig, ctx.fbz.nmtx(iq_fbz), ...
    ctx.fbz.isrtx(:, iq_fbz), ctx.ekin_fbz(:, iq_fbz), 1, ...
    ctx.fbz.mtx{:, iq_fbz}, ctx.gvec, ctx.sys, iq_fbz);
if ctx.use_gpu
    coulg = gpuArray(coulg);
end

if isfield(prepared, 'igpp') && isfield(prepared, 'valid_indices')
    igpp = prepared.igpp;
    valid_indices = prepared.valid_indices;
else
    [igpp, valid_indices] = pre_exact_static_ch( ...
        ctx.fbz, ctx.gvec, kdata.indrk, iq, ctx.use_gpu);
end

block.ik = ik;
block.iq = iq;
block.iq_fbz = iq_fbz;
block.in = in;
block.ispin = ispin;
block.weight = kdata.neq(iq);
block.wfnk = wfnk;
block.wfnkq = wfnkq;
block.idx = idx;
block.fft = ctx.fft;
occ_kq = get_occ(ctx.options, wfnkq.ikq, ispin);
block.occ_kq = occ_kq(1:ctx.nbands);
block.n_cutoff = ctx.fbz.nmtx_cutoff(iq_fbz);
block.coulg = coulg;
block.coulg_cutoff = coulg(1:block.n_cutoff);
block.eps_inv = ctx.eps_inv_fbz{iq_fbz};
block.screened_w = ctx.screened_fbz{iq_fbz};
block.igpp = igpp;
block.valid_indices = valid_indices;
end
