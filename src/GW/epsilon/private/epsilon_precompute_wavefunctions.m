function [wfnk_all, wfnkq_all, fft_all, idx_all] = ...
    epsilon_precompute_wavefunctions(ctx)
%EPSILON_PRECOMPUTE_WAVEFUNCTIONS Precompute epsilon WFN and FFT indices.

wfnk_all = [];
wfnkq_all = [];
fft_all = [];
idx_all = [];
if ~ctx.precompute_wav
    fprintf('No precomputation of wav to save memory.\n');
    return;
end

fprintf('Precomputing wavefunctions...\n');
wfnk_all = cell(ctx.nq, max(cellfun(@(qdata) qdata.nrq, ctx.qdata)));
wfnkq_all = cell(size(wfnk_all));
idx_all.k = cell(ctx.nq, ctx.gr.nf);
idx_all.q = cell(ctx.nq, 1);
idx_all.kq = cell(ctx.nq, ctx.gr.nf);
fft_all = cell(ctx.nq, 1);

precompute_total = sum(cellfun(@(qdata) qdata.nrq, ctx.qdata));
precompute_count = 0;
for iq = 1:ctx.nq
    qdata = ctx.qdata{iq};
    qq = ctx.pol.qpt(iq, :);
    for ik = 1:qdata.nrq
        precompute_count = precompute_count + 1;
        print_progress(precompute_count, precompute_total, ...
            'Message', 'Precompute WFN', 'Task', 'epsilon_precompute');

        rk = ctx.gr.f(qdata.indrk(ik), :);
        wfnk_all{iq, ik} = genwf(rk, ctx.gr, ctx.gvec, ctx.syms, ...
            ctx.sys, ctx.options, ctx.wfc_cutoff, ctx.nbands, ...
            ctx.use_gpu);
        wfnkq_all{iq, ik} = genwf(rk + qq, ctx.gr, ctx.gvec, ...
            ctx.syms, ctx.sys, ctx.options, ctx.wfc_cutoff, ...
            ctx.nbands, ctx.use_gpu);
        [fft_all, idx_all] = epsilon_prefft( ...
            wfnkq_all{iq, ik}, wfnk_all{iq, ik}, iq, ik, ...
            ctx.pol, fft_all, idx_all, ctx.use_gpu);
    end
end
fprintf('\nPrecomputation completed.\n');
end
