function contribution = epsilon_reduced_evaluate(ctx, block)
%EPSILON_REDUCED_EVALUATE Build reduced-basis epsilon contribution.

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
if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
    npairs = numel(block.valence_bands) * ...
        numel(block.conduction_bands);
    isdf_options.rank = ceil(sqrt(npairs) * ctx.eps.isdf.rank_ratio);
end
space = isdf.build_space(left, right, block.idx.q, ...
    size(block.fft.Nfft1), isdf_options);
solver.method = ctx.eps.isdf.reduced_solver;
solver.froErr = ctx.eps.isdf.cauchy_froErr;
solver.MaxIter = ctx.eps.isdf.cauchy_MaxIter;
solver.freq = ctx.pol.freq / ctx.ryd;
polar = isdf.polarizability( ...
    space, block.ev_occ, block.ev_unocc, solver);
contribution.space = space;
contribution.polar = polar;
end
