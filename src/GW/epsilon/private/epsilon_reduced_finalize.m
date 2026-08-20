function eps = epsilon_reduced_finalize(ctx, eps, acc, iq)
%EPSILON_REDUCED_FINALIZE Finalize reduced-basis epsilon outputs.

for ispin = 1:size(acc.rank, 1)
    for ik = 1:size(acc.rank, 2)
        eps.isdf_reduced_rank{iq, ispin, ik} = acc.rank{ispin, ik};
        eps.isdf_reduced_info{iq, ispin, ik} = acc.info{ispin, ik};
    end
end

coulg = coulG_select(ctx.eps, ctx.pol.nmtx(iq), ...
    ctx.pol.isrtx(:, iq), ctx.ekin(:, iq), 0, ...
    ctx.pol.mtx{:, iq}, ctx.gvec, ctx.sys, iq);
if acc.need_full_inverse
    if iq == 1
        eps.inv = cell(ctx.nq, 1);
    end
    eps.inv{iq} = epsilon_invert_pages( ...
        ctx.pol.nmtx(iq), coulg(:), ctx.fact * acc.chi0);
end
if acc.need_screened_w && ~isfield(eps, 'isdf_screened_w')
    eps.isdf_screened_w = cell(ctx.nq, 1);
end
if acc.need_screened_w && ~isempty(acc.zeta_blocks)
    combined_space.zeta_g = cat(2, acc.zeta_blocks{:});
    combined_polar.coeff = ctx.fact * epsilon_page_blkdiag( ...
        acc.coeff_blocks);
    if isa(combined_space.zeta_g, 'gpuArray') && ~isa(coulg, 'gpuArray')
        coulg = gpuArray(coulg);
    end
    eps.isdf_screened_w{iq} = isdf.screened_w( ...
        combined_space, coulg(:), combined_polar);
    eps.isdf_screened_w{iq} = epsilon_gather_screened_w( ...
        eps.isdf_screened_w{iq});
end
end
