function ctx = epsilon_context(sys, options, syms, eps)
ctx.sys = sys;
ctx.options = options;
ctx.syms = syms;
ctx.eps = eps;
ctx.method = gw_resolve_method(eps.isdf, 'epsilon');
ctx.ryd = 13.6056923;
ctx.nbands = eps.nbnd;
ctx.nspin = sys.nspin;
ctx.nspinor = sys.nspinor;
ctx.nq = sys.nkpts;
ctx.wfc_cutoff = 2 * sys.ecut;
ctx.save_mem = eps.save_mem;
ctx.precompute_wav = eps.precompute_wav;
ctx.use_gpu = eps.use_gpu && exist('gpuDevice', 'file');

sigrid = Ggrid(sys, 4 * sys.ecut);
ctx.gvec = Gvector(sigrid, sys);
ctx.gr = fullbz(options, syms, true);
ctx.fact = 4 / (ctx.gr.nf * sys.vol * ctx.nspin * ctx.nspinor);
ctx.pol.qpt = options.kpts;

if eps.freq_dep == 2 && eps.freq_dep_method == 2
    ctx.pol.nfreq_rel = fix(eps.freq_cutoff / eps.delta_freq) + 1;
    ctx.pol.nfreq = ctx.pol.nfreq_rel + eps.nfreq_imag;
    ctx.pol.freq_grid = 0:eps.delta_freq:eps.freq_cutoff;
    tmp = 0:(eps.nfreq_imag - 1);
    ctx.pol.freq_brd = -2i * ctx.ryd * tmp ./ (tmp - eps.nfreq_imag);
    ctx.pol.freq = [ctx.pol.freq_grid, ctx.pol.freq_brd];
else
    ctx.pol.nfreq = 1;
    ctx.pol.freq = 0;
end

ctx.ekin = zeros(ctx.gvec.ng, ctx.nq);
ctx.qdata = cell(ctx.nq, 1);
for iq = 1:ctx.nq
    qq = ctx.pol.qpt(iq, :);
    [ctx.ekin(:, iq), ctx.pol.isrtx(:, iq)] = sortrx( ...
        qq, ctx.gvec.ng, ctx.gvec.mill, sys);
    ctx.pol.nmtx(:, iq) = gcutoff(ctx.gvec.ng, ctx.ekin(:, iq), ...
        ctx.pol.isrtx(:, iq), eps.cutoff);
    ctx.pol.mtx{:, iq} = ctx.gvec.mill( ...
        ctx.pol.isrtx(1:ctx.pol.nmtx(iq), iq), :);
    box_min = zeros(1, 3);
    box_max = zeros(1, 3);
    [box_min, box_max] = get_gvecs_bounds( ...
        ctx.pol.mtx{:, iq}, box_min, box_max);
    ctx.pol.fftgrid{:, iq} = min( ...
        options.wfn_fftgrid + box_max - box_min, options.fftgrid);
    ctx.qdata{iq} = epsilon_qdata(ctx, iq);
end
end
