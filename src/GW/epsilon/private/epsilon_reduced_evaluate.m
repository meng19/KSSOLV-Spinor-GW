function contribution = epsilon_reduced_evaluate(ctx, block)
%EPSILON_REDUCED_EVALUATE Build reduced-basis epsilon contribution.

progress_work = local_progress_work(block);
left = cell(1, ctx.nspinor);
right = cell(1, ctx.nspinor);
for ispinor = 1:ctx.nspinor
    left{ispinor} = isdf.real_component(block.wfnkq, ...
        block.fft.Nfft1, block.idx.kq, block.ispin, ispinor, ...
        block.valence_bands);
    right{ispinor} = isdf.real_component(block.wfnk, ...
        block.fft.Nfft2, block.idx.k, block.ispin, ispinor, ...
        block.conduction_bands);
end
isdf_options = ctx.eps.isdf;
space = isdf.build_space(left, right, block.idx.q, ...
    size(block.fft.Nfft1), isdf_options);
epsilon_progress(block, progress_work * 0.60, ...
    sprintf('E q%d i%d isdf r%d', ...
    block.iq, block.ik, space.rank));
solver.method = ctx.eps.isdf.reduced_solver;
solver.froErr = ctx.eps.isdf.cauchy_froErr;
solver.MaxIter = ctx.eps.isdf.cauchy_MaxIter;
solver.freq = ctx.pol.freq / ctx.ryd;
solver.progress = local_progress_slice(block, progress_work, 0.60, 0.95);
polar = isdf.polarizability( ...
    space, block.ev_occ, block.ev_unocc, solver);
epsilon_progress(block, progress_work * 0.95, ...
    sprintf('E q%d i%d polar %dp', ...
    block.iq, block.ik, numel(block.valence_bands) * ...
    numel(block.conduction_bands)));
contribution.space = space;
contribution.polar = polar;
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

function work = local_progress_work(block)
if isfield(block, 'progress') && isfield(block.progress, 'block_work')
    work = block.progress.block_work;
else
    work = max(1, numel(block.valence_bands) * ...
        numel(block.conduction_bands));
end
end
