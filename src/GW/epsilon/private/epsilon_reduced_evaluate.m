function contribution = epsilon_reduced_evaluate(ctx, block)
%EPSILON_REDUCED_EVALUATE Build reduced-basis epsilon contribution.

progress_work = gw_block_work(block, max(1, ...
    numel(block.valence_bands) * numel(block.conduction_bands)));
[left, left_cache] = isdf.real_components(block.wfnkq, ...
    block.fft.Nfft1, block.idx.kq, block.ispin, ctx.nspinor, ...
    block.valence_bands, [], ctx.eps.isdf.cache_real_wfn);
[right, right_cache] = isdf.real_components(block.wfnk, ...
    block.fft.Nfft2, block.idx.k, block.ispin, ctx.nspinor, ...
    block.conduction_bands, [], ctx.eps.isdf.cache_real_wfn);
if ctx.eps.isdf.cache_real_wfn
    contribution.real_wfn_k = [block.wfnkq.ikq, block.wfnk.ikq];
    contribution.real_wfn = {left_cache, right_cache};
end
% VC products have their own rank controls (rank_vc/rank_ratio_vc), so a
% sigma-oriented rank does not silently determine the epsilon basis.
isdf_options = isdf.options_for_type(ctx.eps.isdf, 'vc');
isdf_options.progress = @(fraction, stage) gw_block_progress(block, ...
    progress_work * 0.60 * fraction, ...
    sprintf('E q%d i%d ISDF %s', block.iq, block.ik, stage));
gw_timer('start', 'Epsilon ISDF interpolation');
space = isdf.build_space(left, right, block.idx.q, ...
    size(block.fft.Nfft1), isdf_options);
gw_timer('stop', 'Epsilon ISDF interpolation');
gw_block_progress(block, progress_work * 0.60, ...
    sprintf('E q%d i%d isdf r%d', ...
    block.iq, block.ik, space.rank));
local_validate_vc_coulomb(ctx, block, space, left, right);
solver.method = ctx.eps.isdf.reduced_solver;
solver.froErr = ctx.eps.isdf.cauchy_froErr;
solver.MaxIter = ctx.eps.isdf.cauchy_MaxIter;
solver.freq = ctx.pol.freq / ctx.ryd;
solver.progress = local_progress_slice(block, progress_work, 0.60, 0.95);
gw_timer('start', 'Epsilon polarizability');
polar = isdf.polarizability( ...
    space, block.ev_occ, block.ev_unocc, solver);
if isfield(space, 'polar_transform')
    polar.coeff = isdf.transform_polar_coeff( ...
        polar.coeff, space.polar_transform);
end

gw_timer('stop', 'Epsilon polarizability');
gw_block_progress(block, progress_work * 0.95, ...
    sprintf('E q%d k%d polar: %d v-c pairs', block.iq, block.ik, ...
    numel(block.valence_bands) * numel(block.conduction_bands)));
contribution.space = space;
contribution.polar = polar;
end

function local_validate_vc_coulomb(ctx, block, space, left, right)
% Epsilon's VC space has no occupied-occupied products, so validate the
% same Coulomb quadratic form as HF exchange on its actual VC products.
if ~ctx.eps.isdf.validate_hf_exchange
    return;
end
weights = repmat(block.occ_vkq(block.valence_bands(:)), 1, ...
    numel(block.conduction_bands));
weights = weights(:);
if ~any(weights)
    return;
end
coulg = coulG_select(ctx.eps, ctx.pol.nmtx(block.iq), ...
    ctx.pol.isrtx(:, block.iq), ctx.ekin(:, block.iq), 1, ...
    ctx.pol.mtx{:, block.iq}, ctx.gvec, ctx.sys, block.iq);
if ctx.use_gpu
    coulg = gpuArray(coulg);
end
validation_options = ctx.eps.isdf;
validation_options.validate_hf_exchange_label = 'VC Coulomb-product';
gw_timer('start', 'Epsilon ISDF VC validation');
report = isdf.validate_hf_exchange(space, left, right, block.idx.q, ...
    size(block.fft.Nfft1), coulg, weights, validation_options);
gw_timer('stop', 'Epsilon ISDF VC validation');
fprintf('ISDF VC interpolation validation for epsilon q%d k%d\n', ...
    block.iq, block.ik);
if isfinite(ctx.eps.isdf.validate_hf_exchange_max_pairs) && ...
        report.pairs_checked < report.pairs_total
    fprintf('ISDF VC validation used the first %d weighted pairs.\n', ...
        report.pairs_checked);
end
end

function progress = local_progress_slice(block, progress_work, start_frac, end_frac)
progress = [];
if ~isfield(block, 'progress') || isempty(block.progress)
    return;
end
progress = block.progress;
progress.completed_before = progress.completed_before + ...
    progress_work * start_frac;
progress.block_work = progress_work * (end_frac - start_frac);
progress.label = sprintf('E q%d i%d s%d', ...
    block.iq, block.ik, block.ispin);
progress.left_bands = block.valence_bands;
progress.right_bands = block.conduction_bands;
end
