function sig = sigma(eps, sig, sys, options, syms)
sig = sigma_set_defaults(sig);
isdf.report_rank('reset');
gw_section_banner('sigma', 'start');
gw_timer('reset');
gw_timer('start', 'Sigma total');
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
precompute_wav = ctx.precompute_wav;

% 添加GPU支持标志
if use_gpu
    fprintf('GPU acceleration enabled for sigma calculation\n');
    gpu_dev = gpuDevice();
    fprintf('Using GPU: %s\n', gpu_dev.Name);
end

ops = sigma_ops(ctx);
if strcmp(ctx.method, 'reduced_basis') || strcmp(ctx.method, 'matrix_elements')
    sigma_isdf_component_cache('reset');
end

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
%%
% Precompute wavefunctions for all k-points and spins
gw_timer('start', 'Sigma wavefunction setup');
[wfnk_all, wfnkq_all, idx_all, igpp, valid_indices] = ...
    sigma_precompute_wavefunctions(ctx);
gw_timer('stop', 'Sigma wavefunction setup');
%%

fprintf('Starting sigma calculation (%d spin channel(s), %d band(s))...\n', ...
    nspin, ndiag);
sigma_q_work = sum(cellfun(@(kdata) kdata.nrk, ctx.kdata));
sigma_block_work = max(1, ctx.nbands);
total_sigma_work = nspin * ndiag * sigma_q_work * sigma_block_work;
current_sigma_work = 0;
sigma_task = 'sigma_main';
progress_percent_step = 10;
progress_update_interval = 5;
print_progress(0, total_sigma_work, ...
    'Message', 'Sigma', ...
    'Task', sigma_task, ...
    'Reset', true, ...
    'StartOnly', true);

for ispin = 1 : nspin
    fprintf('Processing spin %d of %d...\n', ispin, nspin);
    
    for in = ndiag_min : ndiag_max
        for ik = 1 : sig.nkn
            kdata = ctx.kdata{ik};
            rk = kdata.rk;
            
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

                prepared = sigma_prepared_data( ...
                    ctx, ik, iq, wfnk, wfnkq_all, idx_all, ...
                    igpp, valid_indices);
                block = sigma_prepare_block( ...
                    ctx, ik, iq, in, ispin, prepared);
                block.progress = struct( ...
                    'task', sigma_task, ...
                    'completed_before', current_sigma_work, ...
                    'block_work', sigma_block_work, ...
                    'total_work', total_sigma_work, ...
                    'percent_step', progress_percent_step, ...
                    'update_interval', progress_update_interval);
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
                    asx(n_index,ik,ispin) = gw_gather_if_gpu(asxtemp);
                    ax(n_index,ik,ispin) = gw_gather_if_gpu(axtemp);
                    ach(n_index,ik,ispin) = ...
                        0.5 * gw_gather_if_gpu(achtemp);
                    if sig.exact_static_ch
                        achx(n_index,ik,ispin) = ...
                            0.5 * gw_gather_if_gpu(achxtemp);
                    end
                elseif sig.freq_dep == 2
                    asx(n_index,ik,ispin) = ...
                        gw_gather_if_gpu(asxtemp(contribution.iw_lda));
                    asx_freq{n_index,ik,ispin} = ...
                        gw_gather_if_gpu(asxtemp);
                    ax(n_index,ik,ispin) = gw_gather_if_gpu(axtemp);
                    ach(n_index,ik,ispin) = ...
                        gw_gather_if_gpu(achtemp(contribution.iw_lda));
                    ach_freq{n_index,ik,ispin} = ...
                        gw_gather_if_gpu(achtemp);
                    if sig.exact_static_ch
                        achx(n_index,ik,ispin) = ...
                            gw_gather_if_gpu(achxtemp);
                    end
                end
                current_sigma_work = current_sigma_work + ...
                    sigma_block_work;
                print_progress(current_sigma_work, total_sigma_work, ...
                    'Message', sprintf('S b%d k%d q%d complete', ...
                    in, ik, iq), ...
                    'Task', sigma_task, ...
                    'PercentStep', progress_percent_step, ...
                    'UpdateInterval', progress_update_interval);
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

% Retain the energy-resolved terms used by qp.dat writers.  They are in eV
% and have the same diagonal-band layout as eqp0/eqp1.
sig.emf = emf(ndiag_min:ndiag_max, :, :);
sig.ax = real(ax) * ryd;
sig.asx = real(asx) * ryd;
sig.ach = real(ach) * ryd;
sig.achx = real(achx) * ryd;

fprintf('\nCalculation completed.\n');
gw_timer('stop', 'Sigma total');
gw_timer('report', 'Sigma timing information');
gw_section_banner('sigma', 'end');
if strcmp(ctx.method, 'reduced_basis') || strcmp(ctx.method, 'matrix_elements')
    sigma_isdf_component_cache('reset');
end
end
