function sig = sigma(eps, sig, sys, options, syms)
sig = sigma_set_defaults(sig);
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

ndiag = ndiag_max - ndiag_min + 1;
aqsch = cell(nbands, nspin);
asx = zeros([ndiag sys.nkpts nspin]);
ax = zeros([ndiag sys.nkpts nspin]);
ach = zeros([ndiag sys.nkpts nspin]);
achx = zeros([ndiag sys.nkpts nspin]);
sigrid = Ggrid(sys, 4 * sys.ecut);
gvec = Gvector(sigrid,sys);
no_symmetries_q_grid = sig.no_symmetries_q_grid;
sig.qpt = options.kpts;
sig.nkn = sys.nkpts;

% 添加Full frequency支持
if sig.freq_dep == 2 && sig.freq_dep_method == 2
    if sig.freq_grid_shift == 2
        sig.nfreq_grid = 2 * fix(sig.max_freq_eval/sig.delta_freq_eval) + 1; % For Residue of Σ_CH
        sig.freq_grid = 0:sig.delta_freq_eval:2*sig.max_freq_eval;
    end
    sig.nfreq_integral = eps.nfreq; % For integral of Σ_CH
    sig.freq_integral = eps.freq;
    sig.nfreq_integral_imag = eps.nfreq_imag;
    sig.nfreq_integral_real = eps.nfreq - eps.nfreq_imag;
elseif sig.freq_dep == 0
    sig.nfreq_grid = 1;
end
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
    eps_inv_fbz{ik} = eps.inv{gr.indr(ik)}(indt, indt, :);
end

precompute_wav = sig.precompute_wav;
if precompute_wav
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
            wfnk_all{ik} = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff, nbands, use_gpu);
            
            rkq = rk - qq;
            wfnkq_all{iq, ik} = genwf(rkq, gr, gvec, syms, sys, options, wfc_cutoff, nbands, use_gpu);
            
            % 由于FFT格点仅与k, q有关，预计算信息
            idx_all = sigma_prefft(wfnkq_all{iq, ik}, wfnk_all{ik}, fbz.mtx{:, indrk(iq)}, iq, ik, sys, idx_all, use_gpu);
            
            % 如果计算exact_static_ch，由于格点相减仅与k, q有关，预计算信息
            [igpp{iq, ik}, valid_indices{iq, ik}]= pre_exact_static_ch(fbz, gvec, indrk, iq, use_gpu);
        end
    end
else
    fprintf('No precomputation of wav to save memory.\n');
end

%% set up fft_grid
grid_size = [sys.n1, sys.n2, sys.n3];

if use_gpu
    try
        fft.Nfft1 = gpuArray.zeros(grid_size);
        fft.Nfft2 = gpuArray.zeros(grid_size);
        use_gpu = true;
    catch
        warning('GPU memory insufficient for FFT grids. Falling back to CPU.');
        fft.Nfft1 = zeros(grid_size);
        fft.Nfft2 = zeros(grid_size);
        use_gpu = false;
    end
else
    fft.Nfft1 = zeros(grid_size);
    fft.Nfft2 = zeros(grid_size);
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
            
            rk = sig.qpt(ik, :);
            syms_rk = subgrp(rk, syms);
            [nrk, neq, indrk] = irrbz(syms_rk, gr);
            if precompute_wav
                % Use precomputed wfnk
                wfnk = wfnk_all{ik};
            else
                wfnk = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff, nbands, use_gpu);
            end
            if no_symmetries_q_grid
                nrk = gr.nf;
                indrk = (1:nrk);
                neq = ones(1, nrk);
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
                if strcmp(sig.coul_cut, 'spherical_truncation')
                    coulg = coulG_spherical_truncation(fbz.nmtx(1, indrk(iq)), fbz.isrtx(:, indrk(iq)), ekin(:, indrk(iq)), sig.coul_cutoff, 1);
                    
                elseif strcmp(sig.coul_cut, 'cell_box_truncation')
                    if iq > 1
                        error('cell_box_truncation only support one Gamma=0 calculation')
                    end
                    coulg = coulG_cell_box_truncation(fbz.mtx{:, iq}, gvec, sys);
                    
                else
                    error('Unknown truncation schemes for the Coulomb potential: %s. Please choose spherical_truncation or cell_box_truncation.', eps.coul_cut);
                end
                coulg_nocut = fact * coulg;
                coulg_cutoff = coulg_nocut(1 : n_cutoff, 1);
                eps_inv_I = eps_inv - I;
                eps_inv_I_coul = eps_inv_I .* coulg_cutoff';
                if use_gpu
                    eps_inv_I_coul = gpuArray(eps_inv_I_coul);
                    coulg_nocut    = gpuArray(coulg_nocut);
                    coulg_cutoff   = gpuArray(coulg_cutoff);
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
                
                for nn = 1 : nbands
                    aqs{nn, ispin} = getm_sigma(in, nn, wfnkq, wfnk, fft, idx, ispin, nspinor, use_gpu);
                    aqs_nocut = aqs{nn, ispin};
                    aqs_cutoff = aqs{nn, ispin}(1 : n_cutoff, 1);
                    if occ_kq(nn) > 0
                        ax_loc = ax_loc - occ_kq(nn) * sum(abs(aqs_nocut).^2 .* coulg_nocut);
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

fprintf('Calculation completed.\n');
end