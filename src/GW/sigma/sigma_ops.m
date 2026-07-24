function ops = sigma_ops(ctx)
%SIGMA_OPS Select matrix construction while sharing sigma contractions.

ops.name = ctx.method;
switch ctx.method
    case 'direct'
        ops.matrix_elements = @(block) ...
            local_matrix_elements(ctx, block, false);
        ops.contract = @(block, matrix_elements) ...
            local_contract_full(ctx, block, matrix_elements);
    case 'matrix_elements'
        ops.matrix_elements = @(block) ...
            local_matrix_elements(ctx, block, true);
        ops.contract = @(block, matrix_elements) ...
            local_contract_full(ctx, block, matrix_elements);
    case 'reduced_basis'
        error('ISDF:ReducedSigmaNotIntegrated', ...
            'Reduced sigma is integrated in the next migration task.');
end
end

function matrix_elements = local_matrix_elements(ctx, block, use_isdf)
nq = numel(block.idx.q);
if use_isdf
    left = cell(1, ctx.nspinor);
    right = cell(1, ctx.nspinor);
    for ispinor = 1:ctx.nspinor
        left{ispinor} = isdf.real_component(block.wfnk, ...
            block.fft.Nfft1, block.idx.k, block.ispin, ispinor, block.in);
        right{ispinor} = isdf.real_component(block.wfnkq, ...
            block.fft.Nfft2, block.idx.kq, block.ispin, ispinor, ...
            1:ctx.nbands);
    end
    isdf_options = ctx.sig.isdf;
    if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
        isdf_options.rank = ceil( ...
            sqrt(ctx.nbands) * ctx.sig.isdf.rank_ratio);
    end
    gme3 = isdf.matrix_elements(left, right, block.idx.q, ...
        ctx.grid_size, isdf_options);
    matrix_elements = reshape(gme3, nq, ctx.nbands);
    return;
end

if ctx.use_gpu
    matrix_elements = gpuArray.zeros(nq, ctx.nbands);
else
    matrix_elements = zeros(nq, ctx.nbands);
end
for nn = 1:ctx.nbands
    matrix_elements(:, nn) = getm_sigma(block.in, nn, ...
        block.wfnkq, block.wfnk, block.fft, block.idx, block.ispin, ...
        ctx.nspinor, ctx.use_gpu);
end
end

function contribution = local_contract_full(ctx, block, matrix_elements)
eps_inv_I = block.eps_inv - eye(block.n_cutoff);
eps_inv_I_coul = ctx.fact * ...
    (eps_inv_I .* block.coulg_cutoff');
if ctx.use_gpu
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
for nn = 1:ctx.nbands
    aqs_nocut = matrix_elements(:, nn);
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

contribution.asx = asx_loc;
contribution.ax = ax_loc;
contribution.ach = ach_loc;
contribution.achx = achx_loc;
contribution.omega = omega;
contribution.iw_lda = iw_lda;
if ctx.sig.freq_dep == 2
    contribution.asx_freq = asx_loc;
    contribution.ach_freq = ach_loc;
else
    contribution.asx_freq = [];
    contribution.ach_freq = [];
end
contribution.achx_nn = achx_loc_nn;
end
