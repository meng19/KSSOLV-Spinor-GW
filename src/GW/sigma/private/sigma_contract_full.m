function contribution = sigma_contract_full(ctx, block, matrix_elements)
%SIGMA_CONTRACT_FULL Full G-space screened-interaction contraction.

eps_inv_I = block.eps_inv - eye(block.n_cutoff, 'like', block.eps_inv);
eps_inv_I_coul = ctx.fact * ...
    (eps_inv_I .* block.coulg_cutoff');
if ctx.use_gpu && ~isa(eps_inv_I_coul, 'gpuArray')
    eps_inv_I_coul = gpuArray(eps_inv_I_coul);
end

asx_loc = 0;
ax_loc = 0;
ach_loc = 0;
if ctx.sig.freq_dep == 2
    if ctx.use_gpu
        achx_loc_nn = gpuArray.zeros(block.in, ctx.nbands);
    else
        achx_loc_nn = zeros(block.in, ctx.nbands);
    end
else
    achx_loc_nn = [];
end
omega = [];
iw_lda = [];
progress_work = local_progress_work(block);
for nn = 1:ctx.nbands
    aqs_nocut = matrix_elements.gme(:, nn);
    aqs_cutoff = aqs_nocut(1:block.n_cutoff, 1);
    if block.occ_kq(nn) > 0
        ax_loc = ax_loc - block.occ_kq(nn) * ctx.fact * ...
            sum(abs(aqs_nocut).^2 .* block.coulg);
    end
    if ctx.sig.freq_dep == 0
        [asx_loc, ach_loc] = sigma_cohsex(asx_loc, ach_loc, ...
            block.occ_kq(nn), aqs_cutoff, aqs_cutoff, ...
            eps_inv_I_coul);
    elseif ctx.sig.freq_dep == 2
        [asx_loc, ach_loc, achx_loc_nn(block.in, nn), ...
            omega, iw_lda] = sigma_fullfreq(asx_loc, ach_loc, ...
            block.in, nn, block.wfnk.ikq, block.wfnkq.ikq, ...
            block.occ_kq(nn), ctx.options.ev, block.ispin, ...
            aqs_cutoff, aqs_cutoff, eps_inv_I_coul, ctx.sig);
    end
    sigma_progress(block, progress_work * (0.5 + 0.5 * nn / ctx.nbands), ...
        sprintf('S b%d i%d q%d n%d/%d', ...
        block.in, block.ik, block.iq, nn, ctx.nbands));
end

achx_loc = 0;
if ctx.sig.exact_static_ch
    kdata = ctx.kdata{block.ik};
    exact_ch = sigma_cohsex_exact_ch(block.in, block.ispin, ...
        ctx.fbz, kdata.indrk, block.iq, block.aqsch, ...
        eps_inv_I_coul, ctx.sig, block.igpp, block.valid_indices);
    if ctx.sig.freq_dep == 0
        achx_loc = sum(exact_ch, 'all');
    elseif ctx.sig.freq_dep == 2
        achx_loc_nn(block.in, 1) = achx_loc_nn(block.in, 1) + ...
            0.5 * 0.5 * sum(exact_ch, 'all');
        achx_loc = sum(achx_loc_nn(block.in, :), 'all');
    end
end

contribution = sigma_make_contribution( ...
    ctx, asx_loc, ax_loc, ach_loc, achx_loc, ...
    omega, iw_lda, achx_loc_nn);
end

function work = local_progress_work(block)
if isfield(block, 'progress') && isfield(block.progress, 'block_work')
    work = block.progress.block_work;
else
    work = 1;
end
end
