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
nspinor = ctx.nspinor;
use_gpu = ctx.use_gpu;
gvec = ctx.gvec;
gr = ctx.gr;
fact = ctx.fact;
fft = ctx.fft;
ekin = ctx.ekin_fbz;
fbz = ctx.fbz;
eps_inv_fbz = ctx.eps_inv_fbz;
screened_fbz = ctx.screened_fbz; %#ok<NASGU>
precompute_wav = sig.precompute_wav;
no_symmetries_q_grid = sig.no_symmetries_q_grid;

% 添加GPU支持标志
if use_gpu
    fprintf('GPU acceleration enabled for sigma calculation\n');
    gpu_dev = gpuDevice();
    fprintf('Using GPU: %s\n', gpu_dev.Name);
end

if strcmp(ctx.method, 'reduced_basis')
    sig = isdf_sigma_reduced_basis(eps, sig, sys, options, syms);
    return;
end

use_isdf_matrix_elements = strcmp(ctx.method, 'matrix_elements');

ndiag = ndiag_max - ndiag_min + 1;
aqsch = cell(nbands, nspin);
asx = zeros([ndiag sys.nkpts nspin]);
ax = zeros([ndiag sys.nkpts nspin]);
ach = zeros([ndiag sys.nkpts nspin]);
achx = zeros([ndiag sys.nkpts nspin]);
% 添加Full frequency支持
if sig.freq_dep == 2 && sig.freq_dep_method == 2
    % Frequency grids and integration fields are initialized by sigma_context.
elseif sig.freq_dep == 0
    % Static COHSEX uses the single zero-frequency grid from sigma_context.
end
if precompute_wav
    % Precompute wavefunctions for all k-points and spins
    fprintf('Precomputing wavefunctions...\n');
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
            % 预计算neq的GPU版本
            if use_gpu
                neq = gpuArray(neq);
            end
            
            for iq = 1 : nrk
                n_cutoff = fbz.nmtx_cutoff(1, indrk(iq));
                qq = gr.f(indrk(iq), :);
                eps_inv = eps_inv_fbz{indrk(iq)};
                if ~no_symmetries_q_grid
                    [nstar, indst, rqs] = rqstar(syms_rk, qq);
                    if (nstar ~= neq(iq))
                        error('nstar of kpoint %d mismatch', qq)
                    end
                end
                
                %%
                I = eye(fbz.nmtx_cutoff(indrk(iq)));
                coulg = coulG_select(sig, fbz.nmtx(1, indrk(iq)), ...
                    fbz.isrtx(:, indrk(iq)), ekin(:, indrk(iq)), ...
                    1, fbz.mtx{:, indrk(iq)}, gvec, sys, indrk(iq));
                coulg_cutoff = coulg(1 : n_cutoff, 1);
                eps_inv_I = eps_inv - I;
                eps_inv_I_coul = fact * (eps_inv_I .* coulg_cutoff');
                if use_gpu
                    eps_inv_I_coul = gpuArray(eps_inv_I_coul);
                    coulg = gpuArray(coulg);
                    coulg_cutoff = gpuArray(coulg_cutoff);
                end
                
                if precompute_wav
                    %% get wavefunction of k-q from precomputed data
                    wfnkq = wfnkq_all{iq, ik};
                else
                    rkq = rk - qq;
                    wfnkq = genwf(rkq, gr, gvec, syms, sys, options, wfc_cutoff, nbands, use_gpu);
                end
                %% Sum over band nn
                occ_kq = get_occ(options, wfnkq.ikq, ispin);
                if precompute_wav
                    idx.k  = idx_all.k{ik};
                    idx.q  = idx_all.q{iq, ik};
                    idx.kq = idx_all.kq{iq, ik};
                    igpp_tmp = igpp{iq, ik};
                    valid_indices_tmp = valid_indices{iq, ik};
                else
                    idx = sigma_prefft(wfnkq, wfnk, fbz.mtx{:, indrk(iq)}, iq, ik, sys, [], use_gpu);
                    [igpp_tmp, valid_indices_tmp]= pre_exact_static_ch(fbz, gvec, indrk, iq, use_gpu);
                end
                
                asx_loc = 0;
                ax_loc  = 0;
                ach_loc = 0;
                aqs = cell(nbands, nspin);

                if use_isdf_matrix_elements
                    isdf_options = sig.isdf;
                    if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
                        isdf_options.rank = ceil(sqrt(nbands) * sig.isdf.rank_ratio);
                    end
                    aqs_isdf = isdf_sigma_batch(wfnkq, wfnk, fft, idx, ispin, ...
                        nspinor, in, 1:nbands, isdf_options);
                end
                
                for nn = 1 : nbands
                    if use_isdf_matrix_elements
                        aqs{nn, ispin} = aqs_isdf(:, nn);
                    else
                        aqs{nn, ispin} = getm_sigma(in, nn, wfnkq, wfnk, fft, idx, ispin, nspinor, use_gpu);
                    end
                    aqs_nocut = aqs{nn, ispin};
                    aqs_cutoff = aqs{nn, ispin}(1 : n_cutoff, 1);
                    if occ_kq(nn) > 0
                        ax_loc = ax_loc - occ_kq(nn) * fact * sum(abs(aqs_nocut).^2 .* coulg);
                    end
                    if sig.freq_dep == 0
                        [asx_loc, ach_loc] = sigma_cohsex(asx_loc, ach_loc, occ_kq(nn), aqs_cutoff, aqs_cutoff, eps_inv_I_coul);
                    elseif sig.freq_dep == 2
                        [asx_loc, ach_loc, achx_loc_nn(in, nn), omega, iw_lda] = sigma_fullfreq(asx_loc, ach_loc, in, nn, wfnk.ikq, wfnkq.ikq, occ_kq(nn), options.ev, ispin, aqs_cutoff, aqs_cutoff, eps_inv_I_coul, sig);
                        omega_storage(in, ik, ispin, :) = omega;
                        iw_lda_storage(in, ik, ispin) = iw_lda;
                    end
                end
                
                asxtemp = asxtemp + asx_loc * neq(iq);
                axtemp = axtemp + ax_loc * neq(iq);
                achtemp = achtemp + ach_loc * neq(iq);
                
                %% Calculate CH with exact ch correlation
                if sig.exact_static_ch
                    if (indrk(iq) == 1) % Only for q==0
                        aqsch{in, ispin} = aqs{in, ispin};
                    end
                    achx_loc = sigma_cohsex_exact_ch(in, ispin, fbz, indrk, iq, aqsch, eps_inv_I_coul, sig, igpp_tmp, valid_indices_tmp);
                    if sig.freq_dep == 0
                        achxtemp = achxtemp + sum(achx_loc,"all") * neq(iq);
                    elseif sig.freq_dep == 2
                        achx_loc_nn(in, 1) = achx_loc_nn(in, 1) + 0.5 * 0.5 * sum(achx_loc,"all"); % 额外1/2？
                        achx_loc_nn = achx_loc_nn * neq(iq);
                        achxtemp = achxtemp + sum(achx_loc_nn(in, :),"all");
                    end
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
                    asx(n_index,ik,ispin) = asxtemp(iw_lda);
                    asx_freq{n_index,ik,ispin} = asxtemp;
                    ax(n_index,ik,ispin) = axtemp;
                    ach(n_index,ik,ispin) = achtemp(iw_lda);
                    ach_freq{n_index,ik,ispin} = achtemp;
                    if sig.exact_static_ch
                        achx(n_index,ik,ispin) = achxtemp;
                        achx_nn{n_index,ik,ispin} = achx_loc_nn;
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
