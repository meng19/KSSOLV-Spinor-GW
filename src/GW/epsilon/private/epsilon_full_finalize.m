function eps = epsilon_full_finalize(ctx, eps, acc, iq)
%EPSILON_FULL_FINALIZE Build and invert full epsilon pages.

chi0 = acc.chi0;
if ~ctx.save_mem
    for ispin = 1:size(acc.deferred, 1)
        for ik_fbz = 1:size(acc.deferred, 2)
            entries = acc.deferred{ispin, ik_fbz};
            for ientry = 1:numel(entries)
                entry = entries{ientry};
                chi0 = epsilon_add_state_batch( ...
                    chi0, entry{1}, entry{2}, []);
            end
        end
    end
end

chi0 = chi0 * ctx.fact;
nmtx = ctx.pol.nmtx(iq);
coulg = coulG_select(eps, nmtx, ctx.pol.isrtx(:, iq), ...
    ctx.ekin(:, iq), 0, ctx.pol.mtx{:, iq}, ctx.gvec, ctx.sys, iq);

if iq == 1
    eps.inv = cell(ctx.nq, 1);
end
eps.inv{iq} = epsilon_invert_pages(nmtx, coulg(:), chi0);
end
