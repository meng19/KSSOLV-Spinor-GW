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
gw_timer('start', 'Sigma screened kernel');
kernel = isdf.screened_kernel( ...
    block.screened_w, target_zeta, block.coulg_cutoff);
gw_timer('stop', 'Sigma screened kernel');
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
if ctx.sig.freq_dep == 2
    omega = [];
    iw_lda = [];
end
progress_work = gw_block_work(block, 1);
if ctx.sig.freq_dep == 0
    gw_timer('start', 'Sigma static contraction');
    [asx_loc, ax_loc, ach_loc] = local_static_batch_contract( ...
        ctx, block, matrix_elements, kernel(:, :, 1));
    gw_timer('stop', 'Sigma static contraction');
    gw_block_progress(block, progress_work, ...
        sprintf('S b%d k%d q%d static: all %d sum bands', ...
        block.in, block.ik, block.iq, ctx.nbands));
else
for nn = 1:ctx.nbands
    aqs = matrix_elements.gme(:, nn);
    if block.occ_kq(nn) > 0
        aqs_exchange = sigma_exchange_matrix_element(matrix_elements, nn, aqs);
        ax_loc = ax_loc - block.occ_kq(nn) * ctx.fact * ...
            sum(abs(aqs_exchange).^2 .* block.coulg);
    end

    if isfield(matrix_elements, 'coeff')
        coeff = matrix_elements.coeff(:, nn);
    else
        coeff = matrix_elements.space.product_mu(:, nn);
    end
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
    gw_block_progress(block, progress_work * (0.5 + 0.5 * nn / ctx.nbands), ...
        sprintf('S b%d i%d q%d n%d/%d', ...
        block.in, block.ik, block.iq, nn, ctx.nbands));
end
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

if ctx.sig.freq_dep == 2
    contribution = sigma_make_contribution( ...
        ctx, asx_loc, ax_loc, ach_loc, achx_loc, ...
        omega, iw_lda, achx_loc_nn);
else
    contribution = sigma_make_contribution( ...
        ctx, asx_loc, ax_loc, ach_loc, achx_loc);
end
end

function [asx_loc, ax_loc, ach_loc] = local_static_batch_contract( ...
        ctx, block, matrix_elements, kernel)
% Evaluate all static NN contractions as dense matrix products.  For C
% containing one ISDF coefficient column per summation band, the desired
% values are diag(C.' * kernel * conj(C)).

coeff = local_coefficients(matrix_elements, ctx.nbands);
kernel_coeff = kernel * conj(coeff);
screened_values = ctx.fact * sum(coeff .* kernel_coeff, 1);

occ = reshape(block.occ_kq, 1, []);
if isa(screened_values, 'gpuArray') && ~isa(occ, 'gpuArray')
    occ = gpuArray(occ);
end
asx_loc = -sum(occ .* screened_values);
ax_loc = local_batch_exchange(ctx, block, matrix_elements, occ);
ach_loc = sum(screened_values);
end

function ax_loc = local_batch_exchange(ctx, block, matrix_elements, occ)
if isfield(matrix_elements, 'gme_exchange') && ...
        isstruct(matrix_elements.gme_exchange)
    bands = matrix_elements.gme_exchange.bands;
    values = matrix_elements.gme_exchange.values;
    exchange_values = sum(bsxfun(@times, abs(values).^2, block.coulg), 1);
    ax_loc = -ctx.fact * sum(occ(bands) .* exchange_values);
    return;
end
exchange_values = sum(bsxfun(@times, abs(matrix_elements.gme).^2, ...
    block.coulg), 1);
ax_loc = -ctx.fact * sum(occ .* exchange_values);
end

function coeff = local_coefficients(matrix_elements, nbands)
if isfield(matrix_elements, 'coeff')
    coeff = matrix_elements.coeff;
else
    coeff = matrix_elements.space.product_mu(:, 1:nbands);
end
end
