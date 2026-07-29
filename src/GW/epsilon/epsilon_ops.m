function ops = epsilon_ops(ctx)
ops.name = ctx.method;
switch ctx.method
    case 'direct'
        ops.init = @(iq) local_init_full(ctx, iq);
        ops.evaluate = @(block) local_evaluate_full(ctx, block, false);
        ops.accumulate = @(acc, contribution, block) ...
            local_accumulate_full(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) local_finalize_full(ctx, eps, acc, iq);
    case 'matrix_elements'
        ops.init = @(iq) local_init_full(ctx, iq);
        ops.evaluate = @(block) local_evaluate_full(ctx, block, true);
        ops.accumulate = @(acc, contribution, block) ...
            local_accumulate_full(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) local_finalize_full(ctx, eps, acc, iq);
    case 'reduced_basis'
        output_mode = lower(ctx.eps.isdf.output);
        need_full_inverse = any(strcmp( ...
            output_mode, {'full_inverse', 'both'}));
        need_screened_w = any(strcmp( ...
            output_mode, {'screened_w', 'both'}));
        if ~need_full_inverse && ~need_screened_w
            error('ISDF:ReducedEpsilonOutput', ...
                'Unknown ISDF reduced-basis epsilon output "%s".', ...
                ctx.eps.isdf.output);
        end
        ops.init = @(iq) local_init_reduced( ...
            ctx, iq, need_full_inverse, need_screened_w);
        ops.evaluate = @(block) local_evaluate_reduced(ctx, block);
        ops.accumulate = @(acc, contribution, block) ...
            local_accumulate_reduced(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) ...
            local_finalize_reduced(ctx, eps, acc, iq);
end
end

% ---- Reduced-basis epsilon handlers ----

function acc = local_init_reduced( ...
    ctx, iq, need_full_inverse, need_screened_w)
acc.need_full_inverse = need_full_inverse;
acc.need_screened_w = need_screened_w;
if acc.need_full_inverse
    acc.chi0 = zeros(ctx.pol.nmtx(iq), ctx.pol.nmtx(iq), ...
        ctx.pol.nfreq);
else
    acc.chi0 = [];
end
acc.zeta_blocks = {};
acc.coeff_blocks = {};
acc.rank = cell(ctx.nspin, ctx.qdata{iq}.nrq);
acc.info = cell(ctx.nspin, ctx.qdata{iq}.nrq);
end

function contribution = local_evaluate_reduced(ctx, block)
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

function acc = local_accumulate_reduced(~, acc, contribution, block)
coeff_chi = conj(contribution.polar.coeff);
if acc.need_full_inverse
    acc.chi0 = local_accumulate_mapped_chi( ...
        acc.chi0, contribution.space.zeta_g, coeff_chi, block.g_maps);
end
if acc.need_screened_w
    zeta_chi = local_mapped_zeta_chi(contribution.space.zeta_g, ...
        block.g_maps);
    acc.zeta_blocks{end + 1} = zeta_chi;
    acc.coeff_blocks{end + 1} = local_repeat_page_blkdiag( ...
        coeff_chi, numel(block.g_maps));
end
acc.rank{block.ispin, block.ik} = contribution.space.rank;
acc.info{block.ispin, block.ik} = contribution.polar.info;
end

function eps = local_finalize_reduced(ctx, eps, acc, iq)
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
    eps.inv{iq} = local_invert_epsilon_pages( ...
        ctx.pol.nmtx(iq), coulg(:), ctx.fact * acc.chi0);
end
if acc.need_screened_w && ~isfield(eps, 'isdf_screened_w')
    eps.isdf_screened_w = cell(ctx.nq, 1);
end
if acc.need_screened_w && ~isempty(acc.zeta_blocks)
    combined_space.zeta_g = cat(2, acc.zeta_blocks{:});
    combined_polar.coeff = ctx.fact * local_page_blkdiag( ...
        acc.coeff_blocks);
    eps.isdf_screened_w{iq} = isdf.screened_w( ...
        combined_space, coulg(:), combined_polar);
end
end

% ---- Full matrix-element epsilon handlers ----

function acc = local_init_full(ctx, iq)
nmtx = ctx.pol.nmtx(iq);
if ctx.use_gpu
    acc.chi0 = gpuArray.zeros(nmtx, nmtx, ctx.pol.nfreq);
else
    acc.chi0 = zeros(nmtx, nmtx, ctx.pol.nfreq);
end
acc.deferred = cell(ctx.nspin, ctx.gr.nf);
end

function contribution = local_evaluate_full(ctx, block, use_isdf)
if ctx.save_mem && ~use_isdf
    contribution.stream_direct = true;
    return;
end

gme_size = [numel(block.idx.q), numel(block.valence_bands), ...
    numel(block.conduction_bands)];
if ctx.use_gpu
    contribution.gme = gpuArray.zeros(gme_size);
else
    contribution.gme = zeros(gme_size);
end

if use_isdf
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
        isdf_options.rank = ceil(sqrt(numel(block.valence_bands) * ...
            numel(block.conduction_bands)) * ctx.eps.isdf.rank_ratio);
    end
    contribution.gme = isdf.matrix_elements(left, right, ...
        block.idx.q, size(block.fft.Nfft1), isdf_options);
else
    for iv_local = 1:numel(block.valence_bands)
        iv = block.valence_bands(iv_local);
        for ic_local = 1:numel(block.conduction_bands)
            ic = block.conduction_bands(ic_local);
            contribution.gme(:, iv_local, ic_local) = getm_epsilon( ...
                iv, ic, block.wfnkq, block.wfnk, block.fft, block.idx, ...
                block.ispin, ctx.nspinor, ctx.use_gpu);
        end
    end
end

contribution.eden = zeros(numel(block.valence_bands), ...
    numel(block.conduction_bands), ctx.pol.nfreq);
for iv_local = 1:numel(block.valence_bands)
    iv = block.valence_bands(iv_local);
    for ic_local = 1:numel(block.conduction_bands)
        ic = block.conduction_bands(ic_local);
        for ifreq = 1:ctx.pol.nfreq
            freq = ctx.pol.freq(ifreq) / ctx.ryd;
            contribution.eden(iv_local, ic_local, ifreq) = get_eden( ...
                iv, ic, block.wfnkq, block.wfnk, block.ispin, ...
                ctx.options, freq, ctx.eps);
        end
    end
end
end

function acc = local_accumulate_full(ctx, acc, contribution, block)
if isfield(contribution, 'stream_direct')
    acc = local_accumulate_direct_stream(ctx, acc, block);
    return;
end

if ctx.save_mem
    acc.chi0 = local_add_mapped_states( ...
        acc.chi0, contribution.gme, contribution.eden, block.g_maps);
    return;
end

gme = local_mapped_gme(contribution.gme, block.g_maps);
eden = local_repeat_eden(contribution.eden, numel(block.g_maps));
if isempty(acc.deferred{block.ispin, block.ik_fbz})
    acc.deferred{block.ispin, block.ik_fbz} = {{gme, eden}};
else
    acc.deferred{block.ispin, block.ik_fbz}{end + 1} = {gme, eden};
end
end

function acc = local_accumulate_direct_stream(ctx, acc, block)
for iv_local = 1:numel(block.valence_bands)
    iv = block.valence_bands(iv_local);
    for ic_local = 1:numel(block.conduction_bands)
        ic = block.conduction_bands(ic_local);
        vector = getm_epsilon(iv, ic, block.wfnkq, block.wfnk, ...
            block.fft, block.idx, block.ispin, ctx.nspinor, ctx.use_gpu);
        for ifreq = 1:ctx.pol.nfreq
            freq = ctx.pol.freq(ifreq) / ctx.ryd;
            eden = get_eden(iv, ic, block.wfnkq, block.wfnk, ...
                block.ispin, ctx.options, freq, ctx.eps);
            acc.chi0(:, :, ifreq) = local_add_mapped_outer( ...
                acc.chi0(:, :, ifreq), vector, eden, block.g_maps);
        end
    end
end
end

function eps = local_finalize_full(ctx, eps, acc, iq)
chi0 = acc.chi0;
if ~ctx.save_mem
    for ispin = 1:size(acc.deferred, 1)
        for ik_fbz = 1:size(acc.deferred, 2)
            entries = acc.deferred{ispin, ik_fbz};
            for ientry = 1:numel(entries)
                entry = entries{ientry};
                chi0 = local_add_deferred_states(chi0, entry{1}, entry{2});
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

if ctx.use_gpu
    if ~isa(coulg, 'gpuArray')
        coulg_gpu = gpuArray(coulg(:));
    else
        coulg_gpu = coulg(:);
    end
    identity = eye(nmtx, 'gpuArray');
    epsilon_pages = repmat(identity, 1, 1, ctx.pol.nfreq);
    epsilon_pages = epsilon_pages - bsxfun(@times, coulg_gpu, chi0);
    eps.inv{iq} = gather(pagefun(@inv, epsilon_pages));
else
    identity = repmat(eye(nmtx), 1, 1, ctx.pol.nfreq);
    epsilon_pages = identity - bsxfun(@times, coulg(:), chi0);
    inverse_pages = cell(1, size(epsilon_pages, 3));
    for ifreq = 1:size(epsilon_pages, 3)
        inverse_pages{ifreq} = inv(epsilon_pages(:, :, ifreq));
    end
    eps.inv{iq} = cat(3, inverse_pages{:});
end
end

% ---- Full matrix-element accumulation helpers ----

% Accumulate one full matrix-element block in save_mem mode.  The
% irreducible-k star is represented by G-index maps, so each mapped
% contribution is a permutation of the same base chi0 page.
function chi0 = local_add_mapped_states(chi0, gme, eden, g_maps)
gme = reshape(gme, size(gme, 1), []);
eden = reshape(eden, [], size(eden, 3));
for ifreq = 1:size(eden, 2)
    weight = eden(:, ifreq).';
    if isa(gme, 'gpuArray') && ~isa(weight, 'gpuArray')
        weight = gpuArray(weight);
    end
    weighted_gme = bsxfun(@times, gme, weight);
    base_chi = conj(gme) * weighted_gme.';
    for it = 1:numel(g_maps)
        gmap = g_maps{it};
        chi0(:, :, ifreq) = chi0(:, :, ifreq) + ...
            base_chi(gmap, gmap);
    end
end
end

% Streaming direct path: form one outer product for a band pair/frequency,
% then add all equivalent G-reordered pages without recomputing the outer
% product for each star member.
function chi0_page = local_add_mapped_outer( ...
    chi0_page, vector, eden, g_maps)
base_outer = conj(vector) * vector.' * eden;
for it = 1:numel(g_maps)
    gmap = g_maps{it};
    chi0_page = chi0_page + base_outer(gmap, gmap);
end
end

% Deferred full path: consume a block whose mapped star members have already
% been concatenated into the band-pair dimension, allowing one BLAS-style
% matrix product per frequency page in finalize.
function chi0 = local_add_deferred_states(chi0, gme, eden)
gme = reshape(gme, size(gme, 1), []);
eden = reshape(eden, [], size(eden, 3));
for ifreq = 1:size(eden, 2)
    weight = eden(:, ifreq).';
    if isa(gme, 'gpuArray') && ~isa(weight, 'gpuArray')
        weight = gpuArray(weight);
    end
    weighted_gme = bsxfun(@times, gme, weight);
    chi0(:, :, ifreq) = chi0(:, :, ifreq) + ...
        conj(gme) * weighted_gme.';
end
end

% Expand a full gme block over all equivalent G-space reorderings.  The
% mapped blocks are concatenated along the band-pair dimension so deferred
% accumulation can process them as one larger batch.
function mapped = local_mapped_gme(gme, g_maps)
mapped_cells = cell(1, numel(g_maps));
for it = 1:numel(g_maps)
    mapped_cells{it} = gme(g_maps{it}, :, :);
end
mapped = cat(2, mapped_cells{:});
end

% Match local_mapped_gme column order: cat(2, gme maps) becomes repeated
% valence-band slices after reshape, so eden is repeated along its first
% dimension before deferred accumulation flattens it.
function repeated = local_repeat_eden(eden, repeat_count)
if repeat_count == 1
    repeated = eden;
    return;
end
repeated = repmat(eden, repeat_count, 1, 1);
end

% ---- Frequency-page and reduced-basis assembly helpers ----

% Invert each dynamic epsilon page independently.  Static calculations are
% represented as a single page and use the same path.
function inverse_pages = local_invert_epsilon_pages(nmtx, coulg, chi0)
identity = repmat(eye(nmtx), 1, 1, size(chi0, 3));
epsilon_pages = identity - bsxfun(@times, coulg(:), chi0);
inverse_cells = cell(1, size(epsilon_pages, 3));
for ifreq = 1:size(epsilon_pages, 3)
    inverse_cells{ifreq} = inv(epsilon_pages(:, :, ifreq));
end
inverse_pages = cat(3, inverse_cells{:});
end

% Page-wise block diagonal assembly for reduced-basis coefficient blocks.
% Each frequency page gets its own blkdiag because dynamic polarizability
% stores coeff(:,:,ifreq).
function combined = local_page_blkdiag(blocks)
nfreq = size(blocks{1}, 3);
page_cells = cell(1, nfreq);
for ifreq = 1:nfreq
    page_blocks = cell(1, numel(blocks));
    for iblock = 1:numel(blocks)
        page_blocks{iblock} = blocks{iblock}(:, :, ifreq);
    end
    page_cells{ifreq} = blkdiag(page_blocks{:});
end
combined = cat(3, page_cells{:});
end

% Reduced-basis chi0 accumulation.  Build the base reduced chi0 once per
% frequency, then add equivalent star members by G-index permutation.
function chi0 = local_accumulate_mapped_chi(chi0, zeta_g, coeff_chi, g_maps)
base_zeta = conj(zeta_g);
for ifreq = 1:size(coeff_chi, 3)
    base_chi = base_zeta * coeff_chi(:, :, ifreq) * base_zeta';
    for it = 1:numel(g_maps)
        gmap = g_maps{it};
        chi0(:, :, ifreq) = chi0(:, :, ifreq) + base_chi(gmap, gmap);
    end
end
end

% Build the combined zeta projector needed for reduced screened-W output.
% Unlike chi0 accumulation, screened-W must keep the mapped projectors
% explicitly because sigma later contracts in the reduced basis.
function zeta_chi = local_mapped_zeta_chi(zeta_g, g_maps)
zeta_cells = cell(1, numel(g_maps));
for it = 1:numel(g_maps)
    zeta_cells{it} = conj(zeta_g(g_maps{it}, :));
end
zeta_chi = cat(2, zeta_cells{:});
end

% Repeat one reduced coefficient page block for all equivalent star members,
% preserving the page dimension for full-frequency screened-W kernels.
function repeated = local_repeat_page_blkdiag(coeff, repeat_count)
if repeat_count == 1
    repeated = coeff;
    return;
end
blocks = repmat({coeff}, 1, repeat_count);
repeated = local_page_blkdiag(blocks);
end
