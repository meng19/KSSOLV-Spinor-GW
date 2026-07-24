function sig = sigma(eps, sig, sys, options, syms)
sig = sigma_set_defaults(sig);
ctx = sigma_context(eps, sig, sys, options, syms);
sig = ctx.sig;
ryd = ctx.ryd;
nbands = ctx.nbands;
ndiag_min = ctx.band_range(1);
ndiag_max = ctx.band_range(end);
wfc_cutoff = ctx.wfc_cutoff;
nspin = ctx.nspin;
use_gpu = ctx.use_gpu;
gvec = ctx.gvec;
gr = ctx.gr;
fbz = ctx.fbz;
precompute_wav = sig.precompute_wav;

% 添加GPU支持标志
if use_gpu
    fprintf('GPU acceleration enabled for sigma calculation\n');
    gpu_dev = gpuDevice();
    fprintf('Using GPU: %s\n', gpu_dev.Name);
end

ops = sigma_ops(ctx);

ndiag = ndiag_max - ndiag_min + 1;
aqsch = cell(nbands, nspin);
asx = zeros([ndiag sys.nkpts nspin]);
ax = zeros([ndiag sys.nkpts nspin]);
ach = zeros([ndiag sys.nkpts nspin]);
achx = zeros([ndiag sys.nkpts nspin]);
% 添加Full frequency支持
if sig.freq_dep == 2 && sig.freq_dep_method == 2
    % Frequency grids and integration fields are initialized by sigma_context.
    omega_storage = zeros(nbands, sig.nkn, nspin, sig.nfreq_grid);
    iw_lda_storage = zeros(nbands, sig.nkn, nspin);
    asx_freq = cell(ndiag, sig.nkn, nspin);
    ach_freq = cell(ndiag, sig.nkn, nspin);
elseif sig.freq_dep == 0
    % Static COHSEX uses the single zero-frequency grid from sigma_context.
end
if precompute_wav
    % Precompute wavefunctions for all k-points and spins
    fprintf('Precomputing wavefunctions...\n');
    wfnk_all = cell(sig.nkn, 1);
    wfnkq_all = cell(gr.nf, sig.nkn);
    igpp = cell(gr.nf, sig.nkn);
    valid_indices = cell(gr.nf, sig.nkn);
    idx_all.k = cell(sig.nkn, 1);
    idx_all.q = cell(gr.nf, sig.nkn);
    idx_all.kq = cell(gr.nf, sig.nkn); % Dimensions: [iq, ik]
    
    % 计算预计算总数
    precompute_total = 0;
    for ik = 1:sig.nkn
        kdata = ctx.kdata{ik};
        precompute_total = precompute_total + kdata.nrk;
    end
    
    precompute_count = 0;
    for ik = 1:sig.nkn
        kdata = ctx.kdata{ik};
        rk = kdata.rk;
        nrk = kdata.nrk;
        indrk = kdata.indrk;
        for iq = 1:nrk
            qq = gr.f(indrk(iq), :);
            precompute_count = precompute_count + 1;
            print_progress(precompute_count, precompute_total, 'Message', 'Precompute WFN', 'Task', 'sigma_precompute');
            
            wfnk_all{ik} = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff, nbands, use_gpu);
            
            rkq = rk - qq;
            wfnkq_all{iq, ik} = genwf(rkq, gr, gvec, syms, sys, options, wfc_cutoff, nbands, use_gpu);
            
            % 由于FFT格点仅与k, q有关，预计算信息
            idx_all = sigma_prefft(wfnkq_all{iq, ik}, wfnk_all{ik}, fbz.mtx{:, indrk(iq)}, iq, ik, sys, idx_all, use_gpu);
            
            % 如果计算exact_static_ch，由于格点相减仅与k, q有关，预计算信息
            [igpp{iq, ik}, valid_indices{iq, ik}]= pre_exact_static_ch(fbz, gvec, indrk, iq, use_gpu);
        end
    end
    fprintf('\nPrecomputation completed.\n');
else
    fprintf('No precomputation of wav to save memory.\n');
end

fprintf('Starting sigma calculation loop over spins and bands...\n');

for ispin = 1 : nspin
    fprintf('Processing spin %d of %d...\n', ispin, nspin);
    
    for in = ndiag_min : ndiag_max
        fprintf('\n Band %d (index %d/%d)', in, in - ndiag_min + 1, ndiag);
        total_iterations = sig.nkn;
        current_iteration = 0;
        for ik = 1 : sig.nkn
            current_iteration = current_iteration + 1;
            kdata = ctx.kdata{ik};
            rk = kdata.rk;
            
            % 使用print_progress函数更新进度条（每10%刷新）
            print_progress(current_iteration, total_iterations, ...
                'Message', sprintf('Sigma sp:%d ik:%d', ispin, ik), ...
                'Task', sprintf('sigma_spin%d_iband%d', ispin, in), ...
                'PercentStep', 10);
            
            % 在循环开始前初始化GPU变量
            if use_gpu
                asxtemp = gpuArray(0);
                axtemp = gpuArray(0);
                achtemp = gpuArray(0);
                if sig.exact_static_ch
                    achxtemp = gpuArray(0);
                end
            else
                asxtemp = 0;
                axtemp = 0;
                achtemp = 0;
                if sig.exact_static_ch
                    achxtemp = 0;
                end
            end
            
            syms_rk = kdata.syms;
            nrk = kdata.nrk;
            neq = kdata.neq;
            indrk = kdata.indrk;
            
            if precompute_wav
                % Use precomputed wfnk
                wfnk = wfnk_all{ik};
            else
                wfnk = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff, nbands, use_gpu);
            end
            for iq = 1 : nrk
                qq = gr.f(indrk(iq), :);
                if ~sig.no_symmetries_q_grid
                    [nstar, ~, ~] = rqstar(syms_rk, qq);
                    if (nstar ~= neq(iq))
                        if strcmp(ctx.method, 'reduced_basis')
                            error('ISDF:ReducedSigmaStar', ...
                                'Q-point star size does not match its weight.');
                        end
                        error('nstar of kpoint %d mismatch', qq)
                    end
                end

                prepared = struct();
                prepared.wfnk = wfnk;
                if precompute_wav
                    prepared.wfnkq = wfnkq_all{iq, ik};
                    prepared.idx.k = idx_all.k{ik};
                    prepared.idx.q = idx_all.q{iq, ik};
                    prepared.idx.kq = idx_all.kq{iq, ik};
                    prepared.igpp = igpp{iq, ik};
                    prepared.valid_indices = valid_indices{iq, ik};
                end
                block = sigma_prepare_block( ...
                    ctx, ik, iq, in, ispin, prepared);
                matrix_elements = ops.matrix_elements(block);
                if sig.exact_static_ch && block.iq_fbz == 1
                    aqsch{in, ispin} = matrix_elements.gme(:, in);
                end
                block.aqsch = aqsch;
                contribution = ops.contract(block, matrix_elements);
                asxtemp = asxtemp + contribution.asx * block.weight;
                axtemp = axtemp + contribution.ax * block.weight;
                achtemp = achtemp + contribution.ach * block.weight;
                if sig.exact_static_ch
                    achxtemp = achxtemp + contribution.achx * block.weight;
                end
                if sig.freq_dep == 2
                    omega_storage(in, ik, ispin, :) = contribution.omega;
                    iw_lda_storage(in, ik, ispin) = contribution.iw_lda;
                end
                
                n_index = in - ndiag_min + 1;
                if sig.freq_dep == 0
                    asx(n_index,ik,ispin) = asxtemp;
                    ax(n_index,ik,ispin) = axtemp;
                    ach(n_index,ik,ispin) = 0.5 * achtemp;
                    if sig.exact_static_ch
                        achx(n_index,ik,ispin) = 0.5 * achxtemp;
                    end
                elseif sig.freq_dep == 2
                    asx(n_index,ik,ispin) = ...
                        asxtemp(contribution.iw_lda);
                    asx_freq{n_index,ik,ispin} = asxtemp;
                    ax(n_index,ik,ispin) = axtemp;
                    ach(n_index,ik,ispin) = ...
                        achtemp(contribution.iw_lda);
                    ach_freq{n_index,ik,ispin} = achtemp;
                    if sig.exact_static_ch
                        achx(n_index,ik,ispin) = achxtemp;
                    end
                end
            end
        end
    end
end
fprintf('Finalizing calculations...\n');

if sig.exact_static_ch
    if sig.freq_dep == 0
        sig.cor = real(asx + achx) * ryd;
        sig.sig = real(asx + ax + achx) * ryd;
    elseif sig.freq_dep == 2
        sig.cor = real(asx + ach + achx) * ryd;
        sig.sig = real(asx + ax + ach + achx) * ryd;
    end
else
    sig.cor = real(asx + ach) * ryd;
    sig.sig = real(asx + ax + ach) * ryd;
end

emf = ryd * options.ev;
sig = quasi_energy(nspin, ndiag_min, ndiag_max, emf, sys.vxc, sig);
if sig.freq_dep == 2
    sig = get_eqp1(nspin, ndiag_min, ndiag_max, emf, omega_storage, iw_lda_storage, asx_freq, ach_freq, achx, ax, sig);
end

fprintf('\nCalculation completed.\n');
end
