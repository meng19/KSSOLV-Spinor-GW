function eps = epsilon(sys, options, syms, eps)
eps = epsilon_set_defaults(eps);
ctx = epsilon_context(sys, options, syms, eps);
nbands = ctx.nbands;
nspin = ctx.nspin;
nspinor = ctx.nspinor;
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

%% Precompute wavefunctions for all k-points and spins
[wfnk_all, wfnkq_all, fft_all, idx_all] = ...
    epsilon_precompute_wavefunctions(ctx);

%% Main loop
fprintf('Starting main epsilon calculation loop...\n');
ops = epsilon_ops(ctx);
epsilon_block_work = max(1, nvbands * ncbands);
total_epsilon_work = epsilon_block_work * nspin * sum(cellfun( ...
    @(qdata) qdata.nrq, ctx.qdata));
current_epsilon_work = 0;
epsilon_task = 'epsilon_main';
print_progress(0, total_epsilon_work, ...
    'Message', 'Epsilon', ...
    'Task', epsilon_task, ...
    'Reset', true, ...
    'StartOnly', true);

for iq = 1:ctx.nq
    qdata = ctx.qdata{iq};
    qq = ctx.pol.qpt(iq, :);
    nmtx_current = ctx.pol.nmtx(iq);

    fprintf('\n[Epsilon] K-point %2d/%2d | K-vector = (%8.4f, %8.4f, %8.4f) | nmtx = %2d', ...
        iq, ctx.nq, qq(1), qq(2), qq(3), nmtx_current);

    acc = ops.init(iq);
    for ispin = 1:nspin
        for ik = 1:qdata.nrq
            prepared = epsilon_prepared_data( ...
                ctx, iq, ik, wfnk_all, wfnkq_all, fft_all, idx_all);
            block = epsilon_prepare_block(ctx, iq, ik, ispin, prepared);
            block.progress = struct( ...
                'task', epsilon_task, ...
                'completed_before', current_epsilon_work, ...
                'block_work', epsilon_block_work, ...
                'total_work', total_epsilon_work, ...
                'percent_step', 1);
            if isempty(block.valence_bands) || isempty(block.conduction_bands)
                current_epsilon_work = current_epsilon_work + ...
                    epsilon_block_work;
                print_progress(current_epsilon_work, total_epsilon_work, ...
                    'Message', sprintf('Eps q:%d ik:%d s:%d', ...
                    iq, ik, ispin), ...
                    'Task', epsilon_task, ...
                    'PercentStep', 1);
                continue;
            end

            contribution = ops.evaluate(block);
            acc = ops.accumulate(acc, contribution, block);

            current_epsilon_work = current_epsilon_work + ...
                epsilon_block_work;
            print_progress(current_epsilon_work, total_epsilon_work, ...
                'Message', sprintf('Eps q:%d ik:%d s:%d', ...
                iq, ik, ispin), ...
                'Task', epsilon_task, ...
                'PercentStep', 1);
        end
    end
    eps = ops.finalize(eps, acc, iq);
end
eps = epsilon_warn_cauchy_fallback(ctx, eps);

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
