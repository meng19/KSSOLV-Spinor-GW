function [wfnk_all, wfnkq_all, idx_all, igpp, valid_indices] = ...
    sigma_precompute_wavefunctions(ctx)
%SIGMA_PRECOMPUTE_WAVEFUNCTIONS Precompute sigma WFN and index data.

wfnk_all = [];
wfnkq_all = [];
idx_all = [];
igpp = [];
valid_indices = [];
if ~ctx.precompute_wav
    fprintf('No precomputation of wav to save memory.\n');
    return;
end

fprintf('Precomputing wavefunctions...\n');
wfnk_all = cell(ctx.sig.nkn, 1);
wfnkq_all = cell(ctx.gr.nf, ctx.sig.nkn);
igpp = cell(ctx.gr.nf, ctx.sig.nkn);
valid_indices = cell(ctx.gr.nf, ctx.sig.nkn);
idx_all.k = cell(ctx.sig.nkn, 1);
idx_all.q = cell(ctx.gr.nf, ctx.sig.nkn);
idx_all.kq = cell(ctx.gr.nf, ctx.sig.nkn);

precompute_total = 0;
for ik = 1:ctx.sig.nkn
    kdata = ctx.kdata{ik};
    precompute_total = precompute_total + kdata.nrk;
end

precompute_count = 0;
for ik = 1:ctx.sig.nkn
    kdata = ctx.kdata{ik};
    rk = kdata.rk;
    indrk = kdata.indrk;
    for iq = 1:kdata.nrk
        qq = ctx.gr.f(indrk(iq), :);
        precompute_count = precompute_count + 1;
        print_progress(precompute_count, precompute_total, ...
            'Message', 'Precompute WFN', 'Task', 'sigma_precompute');

        wfnk_all{ik} = genwf(rk, ctx.gr, ctx.gvec, ctx.syms, ...
            ctx.sys, ctx.options, ctx.wfc_cutoff, ctx.nbands, ...
            ctx.use_gpu);

        rkq = rk - qq;
        wfnkq_all{iq, ik} = genwf(rkq, ctx.gr, ctx.gvec, ...
            ctx.syms, ctx.sys, ctx.options, ctx.wfc_cutoff, ...
            ctx.nbands, ctx.use_gpu);

        idx_all = sigma_prefft( ...
            wfnkq_all{iq, ik}, wfnk_all{ik}, ...
            ctx.fbz.mtx{:, indrk(iq)}, iq, ik, ctx.sys, ...
            idx_all, ctx.use_gpu);

        [igpp{iq, ik}, valid_indices{iq, ik}] = ...
            pre_exact_static_ch( ...
            ctx.fbz, ctx.gvec, indrk, iq, ctx.use_gpu);
    end
end
fprintf('\nPrecomputation completed.\n');
end
