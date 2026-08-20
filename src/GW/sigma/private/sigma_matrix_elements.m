function matrix_elements = sigma_matrix_elements(ctx, block, use_isdf)
%SIGMA_MATRIX_ELEMENTS Build direct or ISDF sigma matrix elements.

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
    if strcmp(ctx.method, 'reduced_basis')
        space = isdf.build_space(left, right, block.idx.q, ...
            ctx.grid_size, isdf_options);
        gme = reshape(space.zeta_g * space.product_mu, ...
            nq, ctx.nbands);
    else
        gme3 = isdf.matrix_elements(left, right, block.idx.q, ...
            ctx.grid_size, isdf_options);
        gme = reshape(gme3, nq, ctx.nbands);
        space = [];
    end
    matrix_elements.gme = gme;
    matrix_elements.space = space;
    return;
end

if ctx.use_gpu
    gme = gpuArray.zeros(nq, ctx.nbands);
else
    gme = zeros(nq, ctx.nbands);
end
for nn = 1:ctx.nbands
    gme(:, nn) = getm_sigma(block.in, nn, ...
        block.wfnkq, block.wfnk, block.fft, block.idx, block.ispin, ...
        ctx.nspinor, ctx.use_gpu);
end
matrix_elements.gme = gme;
matrix_elements.space = [];
end
