function sig = sigma_gpu(eps, sig, sys, options, syms)
%% Initialize some value
ryd = 13.6056923;
nbands = sig.nbnd;
ndiag_min = sig.ndiag_min;
ndiag_max = sig.ndiag_max;
wfc_cutoff = sys.ecut * 2; % Ha --> Ry
nspin = sys.nspin;
nspinor = sys.nspinor;

% 添加GPU支持标志
use_gpu = sig.use_gpu && exist('gpuDevice', 'file'); % 检查GPU是否可用
if use_gpu
    fprintf('GPU acceleration enabled for sigma calculation\n');
    gpu_dev = gpuDevice();
    fprintf('Using GPU: %s\n', gpu_dev.Name);
end

ndiag = ndiag_max - ndiag_min;
aqsch = cell(nbands, nspin);
asx = zeros([ndiag sys.nkpts nspin]);
ax = zeros([ndiag sys.nkpts nspin]);
ach = zeros([ndiag sys.nkpts nspin]);
achx = zeros([ndiag sys.nkpts nspin]);
sigrid = Ggrid(sys, 4 * sys.ecut);
gvec = Gvector(sigrid,sys);
coul_cutoff = sig.coul_cutoff;
no_symmetries_q_grid = sig.no_symmetries_q_grid;
sig.qpt = options.kpts;
sig.nkn = sys.nkpts;

%%
gr = fullbz(options, syms, true);
fact = 1/(gr.nf * sys.vol);
coulfact = 8 * pi * fact;
eps_inv_fbz = cell([gr.nf 1]);

for ik = 1 : sig.nkn
    rk = sig.qpt(ik, :);
    [ekin(:,ik), sig.isrtx(:,ik)] = sortrx(rk, gvec.ng, gvec.mill, sys);
    sig.nmtx(:,ik) = gcutoff(gvec.ng, ekin(:,ik), sig.isrtx(:,ik), eps.cutoff);
    sig.mtx{:, ik} = gvec.mill(sig.isrtx(1:sig.nmtx(ik), ik), :);
end

for ik = 1 : gr.nf
    rk = gr.f(ik, :);
    [ekin(:,ik), fbz.isrtx(:,ik)] = sortrx(rk, gvec.ng, gvec.mill, sys);
    % fbz is for q-points
    fbz.nmtx(:,ik) = gcutoff(gvec.ng, ekin(:,ik), fbz.isrtx(:,ik), wfc_cutoff);
    fbz.mtx{:, ik} = gvec.mill(fbz.isrtx(1:fbz.nmtx(ik), ik), :);
    fbz.nmtx_cutoff(:,ik) = gcutoff(gvec.ng, ekin(:,ik), fbz.isrtx(:,ik), eps.cutoff);
    fbz.mtx_cutoff{:, ik} = gvec.mill(fbz.isrtx(1:fbz.nmtx_cutoff(ik), ik), :);
    
    %% gmap for eps_inv
    itran = gr.itran(ik);
    qk = gr.r(gr.indr(ik),:) * syms.mtrx{itran,:};
    [~, kgq] = krange(qk, 1e-9);
    for i = 1 : gvec.ng
        isorti(sig.isrtx(i, gr.indr(ik)), 1) = i;
    end
    fbz.isorti(:, ik) = isorti;
    indt = gmap(gvec, syms, sig.nmtx(:,gr.indr(ik)), itran, kgq, fbz.isrtx(:,ik) ,isorti, sys);
    eps_inv_fbz{ik} = eps.inv{gr.indr(ik)}(indt, indt);
end

% Precompute wavefunctions for all k-points and spins
fprintf('Precomputing wavefunctions...\n');
idx_all.k = cell(sig.nkn, 1);
idx_all.q = cell(gr.nf, sig.nkn);
idx_all.kq = cell(gr.nf, sig.nkn); % Dimensions: [iq, ik]

for ik = 1:sig.nkn
    rk = sig.qpt(ik, :);
    syms_rk = subgrp(rk, syms);
    [nrk, neq, indrk] = irrbz(syms_rk, gr);
    if no_symmetries_q_grid
        nrk = gr.nf;
        indrk = (1:nrk);
        neq = ones(1, nrk);
    end
    for iq = 1:nrk
        qq = gr.f(indrk(iq), :);
        wfnk_all{ik} = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff, use_gpu);
        
        rkq = rk - qq;
        wfnkq_all{iq, ik} = genwf(rkq, gr, gvec, syms, sys, options, wfc_cutoff, use_gpu);
        
        % 由于FFT格点仅与k, q有关，预计算信息
        idx_all = sigma_prefft(wfnkq_all{iq, ik}, wfnk_all{ik}, fbz.mtx{:, indrk(iq)}, iq, ik, sys, idx_all, use_gpu);
        
        % 如果计算exact_static_ch，由于格点相减仅与k, q有关，预计算信息
        [igpp{iq, ik}, valid_indices{iq, ik}]= pre_exact_static_ch(fbz, gvec, indrk, iq, use_gpu);
    end
end

%% set up fft_grid
grid_size = [sys.n1, sys.n2, sys.n3];

if use_gpu
    try
        fft.Nfft1 = gpuArray.zeros(grid_size, 'single');
        fft.Nfft2 = gpuArray.zeros(grid_size, 'single');
        use_gpu = true;
    catch
        warning('GPU memory insufficient for FFT grids. Falling back to CPU.');
        fft.Nfft1 = zeros(grid_size, 'single');
        fft.Nfft2 = zeros(grid_size, 'single');
        use_gpu = false;
    end
else
    fft.Nfft1 = zeros(grid_size, 'single');
    fft.Nfft2 = zeros(grid_size, 'single');
end
fft.size = prod(grid_size);

fprintf('Starting sigma calculation loop over spins and bands...\n');

for ispin = 1 : nspin
    fprintf('Processing spin %d of %d...\n', ispin, nspin);
    for in = ndiag_min : ndiag_max
        fprintf('Processing band %d ...\n', in);
        for ik = 1 : sig.nkn
            fprintf('Processing k-point %d of %d...\n', ik, sig.nkn);
            
            % 在循环开始前初始化GPU变量
            if use_gpu
                asxtemp_gpu = gpuArray(0);
                axtemp_gpu = gpuArray(0);
                achtemp_gpu = gpuArray(0);
                if sig.exact_static_ch
                    achxtemp_gpu = gpuArray(0);
                end
            else
                asxtemp = 0;
                axtemp = 0;
                achtemp = 0;
                if sig.exact_static_ch
                    achxtemp = 0;
                end
            end
            
            rk = sig.qpt(ik, :);
            syms_rk = subgrp(rk, syms);
            [nrk, neq, indrk] = irrbz(syms_rk, gr);
            % Use precomputed wfnk
            wfnk = wfnk_all{ik};
            if no_symmetries_q_grid
                nrk = gr.nf;
                indrk = (1:nrk);
                neq = ones(1, nrk);
            end
            
            % 预计算neq的GPU版本
            if use_gpu
                neq_gpu = gpuArray(neq);
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
                coulg = getvcoul(fbz.nmtx(1, indrk(iq)), fbz.isrtx(:, indrk(iq)), ekin(:, indrk(iq)), coul_cutoff, 1);
                coulg_nocut = fact * coulg;
                coulg_cutoff = coulg_nocut(1 : n_cutoff, 1);
                eps_inv_I = eps_inv - I;
                eps_inv_I_coul = eps_inv_I .* coulg_cutoff';
                
                %% get wavefunction of k-q from precomputed data
                wfnkq = wfnkq_all{iq, ik};
                %% Sum over band nn
                occ_kq = get_occ(options, wfnkq.ikq, ispin);
                
                idx.k = idx_all.k{ik};
                idx.q = idx_all.q{iq, ik};
                idx.kq = idx_all.kq{iq, ik};
                
                if use_gpu
                    aqs_gpu = cell(nbands, nspin);
                    % ========== GPU完全加速版本 ==========
                    % 将数据转移到GPU
                    eps_inv_I_coul_gpu = gpuArray(eps_inv_I_coul);
                    coulg_nocut_gpu = gpuArray(coulg_nocut);
                    occ_kq_gpu = gpuArray(occ_kq);
                    
                    % 初始化GPU上的累加变量
                    asx_loc_gpu = gpuArray(zeros(size(eps_inv_I_coul)));
                    ax_loc_gpu = gpuArray(zeros(size(coulg_nocut)));
                    ach_loc_gpu = gpuArray(zeros(size(eps_inv_I_coul)));
                    
                    for nn = 1 : nbands
                        aqs_gpu{nn, ispin} = getm_sigma(in, nn, wfnkq, wfnk, fft, idx, ispin, nspinor, use_gpu);
                        aqs_nocut_gpu = aqs_gpu{nn, ispin};
                        aqs_cutoff_gpu = aqs_gpu{nn, ispin}(1 : n_cutoff, 1);
                        
                        % 在GPU上计算所有矩阵运算
                        aqs_product_gpu = aqs_cutoff_gpu .* aqs_cutoff_gpu';
                        aqs_eps_coul_gpu = aqs_product_gpu .* eps_inv_I_coul_gpu;
                        
                        % 计算abs(aqs_nocut).^2 .* coulg_nocut
                        abs_aqs_sq_gpu = abs(aqs_nocut_gpu).^2;
                        ax_temp_gpu = abs_aqs_sq_gpu .* coulg_nocut_gpu;
                        
                        %% Calculate SX and X
                        if occ_kq(nn) > 0
                            % 使用GPU的累加操作
                            asx_loc_gpu = asx_loc_gpu - occ_kq_gpu(nn) * aqs_eps_coul_gpu;
                            ax_loc_gpu = ax_loc_gpu - occ_kq_gpu(nn) * ax_temp_gpu;
                        end
                        
                        %% Calculate CH without exact ch correlation
                        ach_loc_gpu = ach_loc_gpu + aqs_eps_coul_gpu;
                    end
                    
                    % 在GPU上求和并乘以neq
                    asx_sum_gpu = sum(asx_loc_gpu, "all") * neq_gpu(iq);
                    ax_sum_gpu = sum(ax_loc_gpu, "all") * neq_gpu(iq);
                    ach_sum_gpu = sum(ach_loc_gpu, "all") * neq_gpu(iq);
                    
                    % 累加到总和中
                    asxtemp_gpu = asxtemp_gpu + asx_sum_gpu;
                    axtemp_gpu = axtemp_gpu + ax_sum_gpu;
                    achtemp_gpu = achtemp_gpu + ach_sum_gpu;
                    
                else
                    % ========== CPU版本 ==========
                    asx_loc = 0;
                    ax_loc = 0;
                    ach_loc = 0;
                    aqs = cell(nbands, nspin);
                    for nn = 1 : nbands
                        aqs{nn, ispin} = getm_sigma(in, nn, wfnkq, wfnk, fft, idx, ispin, nspinor, use_gpu);
                        aqs_nocut = aqs{nn, ispin};
                        aqs_cutoff = aqs{nn, ispin}(1 : n_cutoff, 1);
                        
                        aqs_product = aqs_cutoff.* aqs_cutoff';
                        aqs_eps_coul = aqs_product .* eps_inv_I_coul;
                        
                        %% Calculate SX and X
                        if occ_kq(nn) > 0
                            asx_loc = asx_loc - occ_kq(nn) * aqs_eps_coul;
                            ax_loc = ax_loc - occ_kq(nn) * abs(aqs_nocut).^2 .* coulg_nocut;
                        end
                        
                        %% Calculate CH without exact ch correlation
                        ach_loc = ach_loc + aqs_eps_coul;
                    end
                    
                    asxtemp = asxtemp + sum(asx_loc,"all") * neq(iq);
                    axtemp = axtemp + sum(ax_loc,"all") * neq(iq);
                    achtemp = achtemp + sum(ach_loc,"all") * neq(iq);
                end
                
                %% Calculate CH with exact ch correlation
                if sig.exact_static_ch
                    ncouls = fbz.nmtx_cutoff(1, indrk(iq));
                    if (indrk(iq) == 1) % Only for q==0
                        if use_gpu
                            aqsch_gpu{in, ispin} = aqs_gpu{in, ispin};
                        else
                            aqsch{in, ispin} = aqs{in, ispin};
                        end
                    end

                    igpp_tmp = igpp{iq, ik};
                    valid_indices_tmp = valid_indices{iq, ik};
                    
                    if use_gpu
                        aqsch_tmp_gpu = gpuArray(zeros(ncouls));
                        aqsch_tmp_gpu(valid_indices_tmp) = aqsch_gpu{in, ispin}(igpp_tmp(valid_indices_tmp));
                        aqsch_tmp_gpu(~valid_indices_tmp) = 0;
                        achx_loc_gpu = aqsch_tmp_gpu .* eps_inv_I_coul_gpu;
                        achxtemp_gpu = achxtemp_gpu + sum(achx_loc_gpu,"all") * neq_gpu(iq);
                    else
                        aqsch_tmp = zeros(ncouls);
                        aqsch_tmp(valid_indices_tmp) = aqsch{in, ispin}(igpp_tmp(valid_indices_tmp));
                        aqsch_tmp(~valid_indices_tmp) = 0;
                        achx_loc = aqsch_tmp .* eps_inv_I_coul;
                        achxtemp = achxtemp + sum(achx_loc,"all") * neq(iq);
                    end
                end
                
                % 从GPU收集结果或使用CPU结果
                if use_gpu
                    asxtemp = gather(asxtemp_gpu);
                    axtemp = gather(axtemp_gpu);
                    achtemp = gather(achtemp_gpu);
                    if sig.exact_static_ch
                        achxtemp = gather(achxtemp_gpu);
                    end
                end
                
                n_index = in - ndiag_min + 1;
                asx(n_index,ik,ispin) = sum(asxtemp,"all");
                ax(n_index,ik,ispin) = sum(axtemp,"all");
                ach(n_index,ik,ispin) = 0.5 * sum(achtemp,"all");
                if sig.exact_static_ch
                    achx(n_index,ik,ispin) = 0.5 * sum(achxtemp,"all");
                end
            end
        end
    end
    
    fprintf('Finalizing calculations...\n');
    
    if sig.exact_static_ch
        sig.cor = real(asx + achx) * ryd;
        sig.sig = real(asx + ax + achx) * ryd;
    else
        sig.cor = real(asx + ach) * ryd;
        sig.sig = real(asx + ax + ach) * ryd;
    end
    
    emf = ryd * options.ev;
    sig = quasi_energy(nspin, ndiag_min, ndiag_max, emf, sys.vxc, sig);
    
    fprintf('Calculation completed.\n');
end