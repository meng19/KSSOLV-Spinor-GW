function eps = epsilon(sys, options, syms, eps)
eps = epsilon_set_defaults(eps);
ctx = epsilon_context(sys, options, syms, eps);
nbands = ctx.nbands;
nspin = ctx.nspin;
nspinor = ctx.nspinor;
precompute_wav = ctx.precompute_wav;
use_gpu = ctx.use_gpu;

nvbands = eps.nv;
ncbands = eps.nc;

% 添加GPU支持
if use_gpu
    fprintf('GPU acceleration enabled for epsilon calculation\n');
    gpu_dev = gpuDevice();
    fprintf('Using GPU: %s\n', gpu_dev.Name);
    fprintf('Available GPU memory: %.2f GB\n', gpu_dev.AvailableMemory/(1024)^3);
end

fprintf('System parameters: nvbands = %d, ncbands = %d, nbands = %d, nspin = %d, nspinor = %d\n', ...
    nvbands, ncbands, nbands, nspin, nspinor);

% 添加Full frequency支持
if eps.freq_dep == 2 && eps.freq_dep_method == 2
    % 初始化频率相关的存储结构
    fprintf('Initializing full frequency-dependent calculation with %d frequencies\n', ...
        ctx.pol.nfreq);
else
    % 非频率依赖静态COHSEX计算：视为只有一个频率0的特例
    fprintf('Initializing static calculation (frequency = 0)\n');
end

wfnk_all = [];
wfnkq_all = [];
fft_all = [];
idx_all = [];

%% Precompute wavefunctions for all k-points and spins
if precompute_wav
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
                ctx.sys, ctx.options, ctx.wfc_cutoff, nbands, use_gpu);
            wfnkq_all{iq, ik} = genwf(rk + qq, ctx.gr, ctx.gvec, ...
                ctx.syms, ctx.sys, ctx.options, ctx.wfc_cutoff, nbands, use_gpu);
            [fft_all, idx_all] = epsilon_prefft( ...
                wfnkq_all{iq, ik}, wfnk_all{iq, ik}, iq, ik, ...
                ctx.pol, fft_all, idx_all, use_gpu);
        end
    end
    fprintf('\nPrecomputation completed.\n');
else
    fprintf('No precomputation of wav to save memory.\n');
end

%% Main loop
fprintf('Starting main epsilon calculation loop...\n');
ops = epsilon_ops(ctx);

for iq = 1:ctx.nq
    qdata = ctx.qdata{iq};
    qq = ctx.pol.qpt(iq, :);
    nmtx_current = ctx.pol.nmtx(iq);

    fprintf('\n[Epsilon] K-point %2d/%2d | K-vector = (%8.4f, %8.4f, %8.4f) | nmtx = %2d', ...
        iq, ctx.nq, qq(1), qq(2), qq(3), nmtx_current);

    total_blocks_for_k = nspin * qdata.nrq;
    current_blocks_for_k = 0;
    acc = ops.init(iq);
    for ispin = 1:nspin
        for ik = 1:qdata.nrq
            prepared = local_epsilon_prepared_data( ...
                ctx, iq, ik, wfnk_all, wfnkq_all, fft_all, idx_all);
            block = epsilon_prepare_block(ctx, iq, ik, ispin, prepared);
            if isempty(block.valence_bands) || isempty(block.conduction_bands)
                current_blocks_for_k = current_blocks_for_k + 1;
                print_progress(current_blocks_for_k, total_blocks_for_k, ...
                    'Message', 'Epsilon', ...
                    'Task', sprintf('epsilon_k%d', iq), ...
                    'PercentStep', 10);
                continue;
            end

            contribution = ops.evaluate(block);
            acc = ops.accumulate(acc, contribution, block);

            current_blocks_for_k = current_blocks_for_k + 1;
            print_progress(current_blocks_for_k, total_blocks_for_k, ...
                'Message', 'Epsilon', ...
                'Task', sprintf('epsilon_k%d', iq), ...
                'PercentStep', 10);
        end
    end
    eps = ops.finalize(eps, acc, iq);
end

% 存储结果
eps.mtx = ctx.pol.mtx;
eps.nmtx = ctx.pol.nmtx;
eps.nfreq = ctx.pol.nfreq;
eps.freq = ctx.pol.freq;

% 最终清理
if use_gpu
    reset(gpuDevice);
end

fprintf('\nCalculation of epsilon completed successfully.\n');
end

function prepared = local_epsilon_prepared_data( ...
    ctx, iq, ik, wfnk_all, wfnkq_all, fft_all, idx_all)
if ~ctx.precompute_wav
    prepared = [];
    return;
end

prepared.wfnk = wfnk_all{iq, ik};
prepared.wfnkq = wfnkq_all{iq, ik};
prepared.fft = fft_all{iq};
prepared.idx.k = idx_all.k{iq, ik};
prepared.idx.q = idx_all.q{iq};
prepared.idx.kq = idx_all.kq{iq, ik};
end
