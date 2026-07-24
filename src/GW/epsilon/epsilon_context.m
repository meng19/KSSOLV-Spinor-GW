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

if strcmp(ctx.method, 'matrix_elements') && ctx.use_gpu
    error('ISDF:EpsilonGPUUnsupported', ...
        'ISDF epsilon currently supports CPU execution only.');
end
if strcmp(ctx.method, 'reduced_basis') && ctx.use_gpu
    error('ISDF:ReducedEpsilonGPU', ...
        'ISDF reduced-basis epsilon currently supports CPU execution only.');
end
if strcmp(ctx.method, 'reduced_basis') && eps.freq_dep ~= 0
    error('ISDF:ReducedEpsilonFrequency', ...
        'ISDF reduced-basis epsilon requires eps.freq_dep = 0.');
end

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
    ctx.qdata{iq} = local_qdata(ctx, iq);
end
end

function qdata = local_qdata(ctx, iq)
qq = ctx.pol.qpt(iq, :);
syms_q = subgrp(qq, ctx.syms);
[nrq, neq, indrk] = irrbz(syms_q, ctx.gr);

isorti = zeros(ctx.gvec.ng, 1);
for ig = 1:ctx.gvec.ng
    isorti(ctx.pol.isrtx(ig, iq)) = ig;
end

qdata.nrq = nrq;
qdata.neq = neq;
qdata.indrk = indrk;
qdata.syms = syms_q;
qdata.g_maps = cell(nrq, 1);
qdata.rqs = cell(nrq, 1);
for ik = 1:nrq
    rk = ctx.gr.f(indrk(ik), :);
    [nstar, indst, rqs] = rqstar(syms_q, rk);
    if nstar ~= neq(ik)
        if strcmp(ctx.method, 'reduced_basis')
            error('ISDF:ReducedEpsilonStar', ...
                'K-point star size does not match its irreducible weight.');
        end
        error('nstar of kpoint %d mismatch', rk);
    end
    qdata.rqs{ik} = rqs;
    qdata.g_maps{ik} = cell(nstar, 1);
    qdata.g_maps{ik}{1} = (1:ctx.pol.nmtx(iq)).';
    for it = 2:nstar
        itran = syms_q.indsub(indst(it));
        kgq = -syms_q.kgzero(indst(it), :);
        qdata.g_maps{ik}{it} = gmap( ...
            ctx.gvec, ctx.syms, ctx.pol.nmtx(iq), itran, kgq, ...
            ctx.pol.isrtx(:, iq), isorti, ctx.sys);
    end
end
end
