function eps = epsilon_reduced_finalize(ctx, eps, acc, iq)
%EPSILON_REDUCED_FINALIZE Finalize reduced-basis epsilon outputs.

for ispin = 1:size(acc.rank, 1)
    for ik = 1:size(acc.rank, 2)
        eps.isdf_reduced_rank{iq, ispin, ik} = acc.rank{ispin, ik};
        eps.isdf_reduced_info{iq, ispin, ik} = acc.info{ispin, ik};
    end
end

if ~isempty(acc.real_wfn)
    if isfield(eps, 'isdf_real_wfn') && ~isempty(eps.isdf_real_wfn)
        for ik = 1:size(acc.real_wfn, 1)
            for ispin = 1:size(acc.real_wfn, 2)
                for ispinor = 1:size(acc.real_wfn, 3)
                    if ~isempty(acc.real_wfn{ik, ispin, ispinor})
                        eps.isdf_real_wfn{ik, ispin, ispinor} = ...
                            acc.real_wfn{ik, ispin, ispinor};
                    end
                end
            end
        end
    else
        eps.isdf_real_wfn = acc.real_wfn;
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
    gw_timer('start', 'Epsilon screened W');
    if numel(acc.zeta_blocks) == 1
        combined_space.zeta_g = acc.zeta_blocks{1};
    else
        combined_space.zeta_g = cat(2, acc.zeta_blocks{:});
    end
    combined_polar.coeff = epsilon_page_blkdiag(acc.coeff_blocks);
    % The combined arrays now own the only references needed by screened W.
    % Release the original blocks before the reduced matrix factorizations.
    acc.zeta_blocks = {};
    acc.coeff_blocks = {};
    if isa(combined_space.zeta_g, 'gpuArray') && ~isa(coulg, 'gpuArray')
        coulg = gpuArray(coulg);
    end
    screened = isdf.screened_w( ...
        combined_space, coulg(:), combined_polar, ...
        ctx.eps.isdf.screened_w_storage);
    % In minimal storage mode these inputs are no longer needed.  Release
    % them before gathering k_mu and zeta_g to the host.
    clear combined_space combined_polar;
    screened = epsilon_gather_screened_w(screened);
    eps.isdf_screened_w{iq} = screened;
    gw_timer('stop', 'Epsilon screened W');
end
end
