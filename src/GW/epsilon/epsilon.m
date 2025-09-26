function eps = epsilon(sys, options, syms, eps)
%% Initialize some value
nvbands = eps.nv;
ncbands = eps.nc;
nbands = eps.nbnd;
nspin = sys.nspin;
nspinor = sys.nspinor;
wfc_cutoff = sys.ecut * 2; % Ha --> Ry

save_mem = eps.save_mem;

% 添加GPU支持
use_gpu = eps.use_gpu && exist('gpuDevice', 'file');
if use_gpu
    fprintf('GPU acceleration enabled for epsilon calculation\n');
    gpu_dev = gpuDevice();
    fprintf('Using GPU: %s\n', gpu_dev.Name);
    fprintf('Available GPU memory: %.2f GB\n', gpu_dev.AvailableMemory/(1024)^3);
end

fprintf('System parameters: nvbands = %d, ncbands = %d, nbands = %d, nspin = %d, nspinor = %d\n', nvbands, ncbands, nbands, nspin, nspinor);

sigrid = Ggrid(sys, 4*sys.ecut);
gvec = Gvector(sigrid, sys);
pol.qpt = options.kpts;

%%
gr = fullbz(options, syms, true);
ekin = zeros(gvec.ng, sys.nkpts);
for iq = 1:sys.nkpts
    qq = pol.qpt(iq,:);
    [ekin(:,iq), pol.isrtx(:,iq)] = sortrx(qq, gvec.ng, gvec.mill, sys);
    pol.nmtx(:,iq) = gcutoff(gvec.ng, ekin(:,iq), pol.isrtx(:,iq), eps.cutoff);
    pol.mtx{:, iq} = gvec.mill(pol.isrtx(1:pol.nmtx(iq), iq), :);
    
    %% get fftsize for calculating M matrix
    eps_box_min = zeros([1 3]);
    eps_box_max = zeros([1 3]);
    [eps_box_min, eps_box_max] = get_gvecs_bounds(pol.mtx{:, iq}, eps_box_min, eps_box_max);
    pol.fftgrid{:, iq} = min((options.wfn_fftgrid + eps_box_max - eps_box_min), options.fftgrid);
end

% 清空不必要的变量
clear sigrid;

eps_tmp = cell(sys.nkpts, 1);
eps_inv = cell(sys.nkpts, 1);
fact = 4 / (gr.nf * sys.vol * nspin * nspinor);

%% Precompute wavefunctions for all k-points and spins
fprintf('Precomputing wavefunctions...\n');
idx_all.k = cell(sys.nkpts, gr.nf);
idx_all.q = cell(sys.nkpts, 1);
idx_all.kq = cell(sys.nkpts, gr.nf); % Dimensions: [iq, ik]
fft_all = cell(sys.nkpts, 1);

for iq = 1:sys.nkpts
    qq = pol.qpt(iq,:);
    syms_qq = subgrp(qq, syms);
    [nrq, neq, indrk] = irrbz(syms_qq, gr);
    for ik = 1:nrq
        rk = gr.f(indrk(ik),:);
        wfnk_all{iq, ik} = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff, use_gpu);
        
        rkq = rk + qq;
        wfnkq_all{iq, ik} = genwf(rkq, gr, gvec, syms, sys, options, wfc_cutoff, use_gpu);
        
        % 由于FFT格点仅与k, q有关，预计算信息
        [fft_all, idx_all] = epsilon_prefft(wfnkq_all{iq, ik}, wfnk_all{iq, ik}, iq, ik, pol, fft_all, idx_all, use_gpu);
    end
end
%% Main loop

for iq = 1:sys.nkpts
    fprintf('\nProcessing k-point %d/%d\n', iq, sys.nkpts);
    
    qq = pol.qpt(iq,:);
    syms_qq = subgrp(qq, syms);
    [nrq, neq, indrk] = irrbz(syms_qq, gr);
    
    nmtx_current = pol.nmtx(iq);
    fprintf('nmtx for current k-point: %d\n', nmtx_current);
    
    % 初始化 chi0 累加器
    if use_gpu
        chi0_sum = gpuArray(zeros(nmtx_current, nmtx_current));
    else
        chi0_sum = zeros(nmtx_current, nmtx_current);
    end
    
    % 预计算所有不可约k点的映射关系
    indt_cell = cell(nrq, 1);
    for ik = 1:nrq
        rk = gr.f(indrk(ik),:);
        [nstar, indst, rqs] = rqstar(syms_qq, rk);
        
        if (nstar ~= neq(ik))
            error('nstar of kpoint %d mismatch', rk);
        end
        
        indt_cell{ik} = cell(nstar, 1);
        for it = 1:nstar
            itran = syms_qq.indsub(indst(it));
            kgq = -syms_qq.kgzero(indst(it),:);
            isorti = zeros(gvec.ng, 1);
            for i = 1:gvec.ng
                isorti(pol.isrtx(i, iq)) = i;
            end
            indt_cell{ik}{it} = gmap(gvec, syms, nmtx_current, itran, kgq, pol.isrtx(:,iq), isorti, sys);
        end
    end
    
    for ispin = 1 : nspin
        for ik = 1 : nrq
            % 初始化临时 chi0
            if use_gpu
                chi0_tmp = gpuArray(zeros(nmtx_current, nmtx_current));
            else
                chi0_tmp = zeros(nmtx_current, nmtx_current);
            end
            fprintf('  Irreducible k-point %d/%d...\n', ik, nrq);
            rk = gr.f(indrk(ik),:);
            [nstar, ~, rqs] = rqstar(syms_qq, rk);
            
            % 读取波函数
            wfnk = wfnk_all{iq, ik};
            wfnkq = wfnkq_all{iq, ik};
            
            % 读取FFT网格
            idx.k = idx_all.k{iq, ik};
            idx.q = idx_all.q{iq};
            idx.kq = idx_all.kq{iq, ik};
            fft = fft_all{iq};
            
            % 获取能带信息
            occ_vkq = get_occ(options, wfnkq.ikq, ispin);
            no_v = sum(occ_vkq > 0);
            occ_ck = get_occ(options, wfnk.ikq, ispin);
            no_c_start = sum(occ_ck > 0) + 1;
            no_c = nbands - no_c_start + 1;
            
            if no_v == 0 || no_c == 0
                continue;
            end
            
            fprintf('    Calculating M matrix elements for ik %d spin %d (v=%d, c=%d)...\n', ik, ispin, no_v, no_c);
            
            % 处理所有价带和导带
            fprintf('开始计算矩阵元...\n');
            for iv = 1:no_v
                % 每10个价带显示一次进度
                if mod(iv, 10) == 1 || iv == no_v
                    fprintf('处理价带 %d/%d', iv, no_v);
                    if no_c > 1
                        fprintf('，导带范围: %d-%d\n', no_c_start, no_c_start + no_c - 1);
                    else
                        fprintf('\n');
                    end
                end
                
                for ic_idx = 1:no_c
                    ic = no_c_start + ic_idx - 1; % 转换为全局能带索引
                    
                    % 计算矩阵元
                    gme_temp = getm_epsilon(iv, ic, wfnkq, wfnk, fft, idx, ispin, nspinor, options, use_gpu);
                    
                    if save_mem
                        % 立即累加到 chi0
                        chi0_tmp = chi0_tmp - conj(gme_temp) * gme_temp.';
                    else
                        % 存储数据（使用局部索引），用于后续使用矩阵乘法加速
                        gme_storage(:, iv, ic_idx, indrk(ik)) = gme_temp;
                    end
                    
                    % 处理简并k点
                    if nstar > 1
                        for it = 2:nstar
                            gme_temp_degen = gme_temp(indt_cell{ik}{it});
                            if save_mem
                                chi0_tmp = chi0_tmp - conj(gme_temp_degen) * gme_temp_degen.';
                            else
                                k_degenerate = rqs(it, :);
                                [~, ik_degenerate] = ismember(k_degenerate, gr.f, 'rows');
                                if ik_degenerate > 0
                                    gme_storage(:, iv, ic_idx, ik_degenerate) = gme_temp_degen;
                                end
                            end
                        end
                    end
                end
            end
            if save_mem
                chi0_sum = chi0_sum + chi0_tmp; % Sum over k and spin
            end
        end
        if ~save_mem
            chi0_sum = chi0_sum + get_chi0(gme_storage); % Sum over k, band and spin
            clear gme_storage
        end
    end
    
    clear chi0_tmp;
    
    % 应用缩放因子
    chi0_sum = chi0_sum * fact;
    
    % 计算Coulomb势
    coulg = getvcoul(nmtx_current, pol.isrtx(:, iq), ekin(:, iq), eps.coul_cutoff, 0);
    
    % 计算epsilon矩阵
    if use_gpu
        coulg_gpu = gpuArray(coulg);
        eps_tmp_gpu = eye(nmtx_current, 'gpuArray');
        eps_tmp_gpu = eps_tmp_gpu - coulg_gpu .* chi0_sum;
        eps_inv_gpu = inv(eps_tmp_gpu);
        
        % 传输回CPU
        eps_tmp{iq} = gather(eps_tmp_gpu);
        eps_inv{iq} = gather(eps_inv_gpu);
        
        clear chi0_sum coulg_gpu eps_tmp_gpu eps_inv_gpu gme_storage;
        wait(gpuDevice);
    else
        eps_tmp{iq} = eye(nmtx_current);
        eps_tmp{iq} = eps_tmp{iq} - coulg .* chi0_sum;
        eps_inv{iq} = inv(eps_tmp{iq});
        
        clear chi0_sum coulg gme_storage;
    end
end

% 存储结果
eps.inv = eps_inv;
eps.mtx = pol.mtx;
eps.nmtx = pol.nmtx;

% 最终清理
if use_gpu
    reset(gpuDevice);
end

fprintf('Calculation completed successfully.\n');
end