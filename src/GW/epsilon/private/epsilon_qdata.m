function qdata = epsilon_qdata(ctx, iq)
%EPSILON_QDATA Build irreducible k data for one epsilon q-point.

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
