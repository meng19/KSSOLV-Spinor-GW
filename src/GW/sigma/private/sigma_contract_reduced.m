function contribution = sigma_contract_reduced(ctx, block, matrix_elements)
%SIGMA_CONTRACT_REDUCED Reduced-basis screened-interaction contraction.

if isempty(block.screened_w)
    if isempty(block.eps_inv)
        error('ISDF:ReducedSigmaMissingQPoint', ...
            'No screened interaction is available for full-BZ q-point %d.', ...
            block.iq_fbz);
    end
    contribution = sigma_contract_full(ctx, block, matrix_elements);
    return;
end

if size(block.screened_w.zeta_g, 1) ~= block.n_cutoff
    error('ISDF:ReducedSigmaScreenedSize', ...
        'Reduced screened interaction does not match sigma cutoff.');
end

target_zeta = matrix_elements.space.zeta_g(1:block.n_cutoff, :);
kernel = isdf.screened_kernel( ...
    block.screened_w, target_zeta, block.coulg_cutoff);
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
    aqs = matrix_elements.gme(:, nn);
    if block.occ_kq(nn) > 0
        ax_loc = ax_loc - block.occ_kq(nn) * ctx.fact * ...
            sum(abs(aqs).^2 .* block.coulg);
    end

    coeff = matrix_elements.space.product_mu(:, nn);
    if ctx.sig.freq_dep == 0
        kernel_static = kernel(:, :, 1);
        screened_value = ctx.fact * isdf.screened_contract( ...
            kernel_static, coeff);
        if block.occ_kq(nn) > 0
            asx_loc = asx_loc - block.occ_kq(nn) * screened_value;
        end
        ach_loc = ach_loc + screened_value;
    elseif ctx.sig.freq_dep == 2
        [asx_loc, ach_loc, achx_loc_nn(block.in, nn), ...
            omega, iw_lda] = sigma_fullfreq(asx_loc, ach_loc, ...
            block.in, nn, block.wfnk.ikq, block.wfnkq.ikq, ...
            block.occ_kq(nn), ctx.options.ev, block.ispin, ...
            coeff, coeff, ctx.fact * kernel, ctx.sig);
    end
    sigma_progress(block, progress_work * (0.5 + 0.5 * nn / ctx.nbands), ...
        sprintf('S b%d i%d q%d n%d/%d', ...
        block.in, block.ik, block.iq, nn, ctx.nbands));
end

achx_loc = 0;
if ctx.sig.exact_static_ch
    screened_matrix = ctx.fact * isdf.screened_kernel( ...
        block.screened_w, [], block.coulg_cutoff);
    kdata = ctx.kdata{block.ik};
    exact_ch = sigma_cohsex_exact_ch(block.in, block.ispin, ...
        ctx.fbz, kdata.indrk, block.iq, block.aqsch, ...
        screened_matrix, ctx.sig, block.igpp, block.valid_indices);
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
