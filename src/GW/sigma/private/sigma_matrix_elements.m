function matrix_elements = sigma_matrix_elements(ctx, block, use_isdf)
%SIGMA_MATRIX_ELEMENTS Build direct or ISDF sigma matrix elements.

nq = numel(block.idx.q);
progress_work = local_progress_work(block);
if use_isdf
    if isfield(ctx.sig.isdf, 'global_nn_space') && ...
            ctx.sig.isdf.global_nn_space
        matrix_elements = local_global_nn_matrix_elements( ...
            ctx, block, progress_work);
    else
        left = local_left_components(ctx, block, block.in);
        right = local_right_components(ctx, block);
        % The generic sigma path needs the NN space for screened and COH
        % terms.  rank_nn is independent of epsilon's VC rank.
        isdf_options = isdf.options_for_type(ctx.sig.isdf, 'nn');
    if strcmp(ctx.method, 'reduced_basis')
            gw_timer('start', 'Sigma ISDF space');
            space = isdf.build_space(left, right, block.idx.q, ...
                ctx.grid_size, isdf_options);
            gw_timer('stop', 'Sigma ISDF space');
            gme = reshape(space.zeta_g * space.product_mu, ...
                nq, ctx.nbands);
        else
            gme3 = isdf.matrix_elements(left, right, block.idx.q, ...
                ctx.grid_size, isdf_options);
            gme = reshape(gme3, nq, ctx.nbands);
            space = [];
        end
        sigma_progress(block, progress_work * 0.5, ...
            sprintf('S b%d i%d q%d me %d/%d', ...
            block.in, block.ik, block.iq, ctx.nbands, ctx.nbands));
        matrix_elements.gme = gme;
        matrix_elements.space = space;
    end
    if local_uses_vn_exchange(ctx)
        matrix_elements.gme_exchange = local_vn_matrix_elements(ctx, block);
    else
        matrix_elements.gme_exchange = [];
    end
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
    sigma_progress(block, progress_work * 0.5 * nn / ctx.nbands, ...
        sprintf('S b%d i%d q%d me %d/%d', ...
        block.in, block.ik, block.iq, nn, ctx.nbands));
end
matrix_elements.gme = gme;
matrix_elements.space = [];
end

function matrix_elements = local_global_nn_matrix_elements( ...
        ctx, block, progress_work)
% Build the NN interpolation space once for all requested diagonal bands.
% This mirrors the vc/vn/nn-space reuse in the reference implementation.
left_bands = ctx.band_range(:).';
key = sprintf('nn-space-k%d-q%d-s%d-n%d-b%d-%d', ...
    block.ik, block.iq, block.ispin, ctx.nbands, ...
    left_bands(1), left_bands(end));
[entry, hit] = sigma_isdf_component_cache('get', key);
if ~hit
    left = local_left_components(ctx, block, left_bands);
    right = local_right_components(ctx, block);
    isdf_options = isdf.options_for_type(ctx.sig.isdf, 'nn');
    if ~strcmp(ctx.method, 'reduced_basis')
        error('ISDF:GlobalNNSpaceRequiresReducedBasis', ...
            'global_nn_space requires the reduced_basis ISDF algorithm.');
    end
    gw_timer('start', 'Sigma ISDF space');
    space = isdf.build_space(left, right, block.idx.q, ...
        ctx.grid_size, isdf_options);
    gw_timer('stop', 'Sigma ISDF space');
    nq = numel(block.idx.q);
    gme_all = reshape(space.zeta_g * space.product_mu, ...
        nq, numel(left_bands), ctx.nbands);
    entry = struct('left_bands', left_bands, 'space', space, ...
        'gme_all', gme_all);
    sigma_isdf_component_cache('put', key, entry);
end

left_index = find(entry.left_bands == block.in, 1);
if isempty(left_index)
    error('ISDF:GlobalNNSpaceBand', ...
        'Requested diagonal band %d is absent from global NN space.', block.in);
end
matrix_elements.gme = reshape(entry.gme_all(:, left_index, :), ...
    numel(block.idx.q), ctx.nbands);
matrix_elements.space = entry.space;
% product_mu is ordered as (left band, right band).  Unlike gme_all, it
% is not reshaped above, so retain the coefficients for this target band
% explicitly for the reduced screened-interaction contraction.
matrix_elements.coeff = entry.space.product_mu(:, ...
    left_index:numel(entry.left_bands):end);
sigma_progress(block, progress_work * 0.5, ...
    sprintf('S b%d i%d q%d global-nn %d/%d', ...
    block.in, block.ik, block.iq, ctx.nbands, ctx.nbands));
end

function gme_exchange = local_vn_matrix_elements(ctx, block)
%VN is the occupied-state/target-state product space used by bare exchange.
occupied_bands = find(block.occ_kq > 0).';
if isempty(occupied_bands)
    gme_exchange = [];
    return;
end

if isfield(ctx.sig.isdf, 'global_vn_space') && ctx.sig.isdf.global_vn_space
    left_bands = ctx.band_range(:).';
else
    left_bands = block.in;
end
key = sprintf('vn-space-k%d-q%d-s%d-b%d-%d-o%s', ...
    block.ik, block.iq, block.ispin, left_bands(1), left_bands(end), ...
    sprintf('%d_', occupied_bands));
[entry, hit] = sigma_isdf_component_cache('get', key);
if ~hit
    left = local_left_components(ctx, block, left_bands);
    right = local_selected_right_components(ctx, block, occupied_bands);
    isdf_options = isdf.options_for_type(ctx.sig.isdf, 'vn');
    nq = numel(block.idx.q);
    if strcmp(ctx.method, 'reduced_basis')
        gw_timer('start', 'Sigma ISDF space');
        space = isdf.build_space(left, right, block.idx.q, ...
            ctx.grid_size, isdf_options);
        gw_timer('stop', 'Sigma ISDF space');
        gme_all = reshape(space.zeta_g * space.product_mu, nq, ...
            numel(left_bands), numel(occupied_bands));
    else
        gme_all = reshape(isdf.matrix_elements(left, right, block.idx.q, ...
            ctx.grid_size, isdf_options), nq, numel(left_bands), ...
            numel(occupied_bands));
    end
    entry = struct('left_bands', left_bands, ...
        'occupied_bands', occupied_bands, 'gme_all', gme_all);
    sigma_isdf_component_cache('put', key, entry);
end

left_index = find(entry.left_bands == block.in, 1);
if isempty(left_index)
    error('ISDF:VNSpaceBand', ...
        'Requested diagonal band %d is absent from VN space.', block.in);
end
gme_exchange = struct('bands', entry.occupied_bands, ...
    'values', reshape(entry.gme_all(:, left_index, :), ...
    numel(block.idx.q), []));
end

function left = local_left_components(ctx, block, bands)
left = local_components(ctx, block.wfnk, block.fft.Nfft1, ...
    block.idx.k, block.ispin, bands);
end

function right = local_right_components(ctx, block)
key = sprintf('k%d-q%d-s%d-n%d', ...
    block.ik, block.iq, block.ispin, ctx.nbands);
[right, hit] = sigma_isdf_component_cache('get', key);
if hit
    return;
end

right = local_components(ctx, block.wfnkq, block.fft.Nfft2, ...
    block.idx.kq, block.ispin, 1:ctx.nbands);
sigma_isdf_component_cache('put', key, right);
end

function right = local_selected_right_components(ctx, block, bands)
key = sprintf('right-k%d-q%d-s%d-b%s', ...
    block.ik, block.iq, block.ispin, sprintf('%d_', bands));
[right, hit] = sigma_isdf_component_cache('get', key);
if hit
    return;
end
right = local_components(ctx, block.wfnkq, block.fft.Nfft2, ...
    block.idx.kq, block.ispin, bands);
sigma_isdf_component_cache('put', key, right);
end

function components = local_components(ctx, wfn, fft_template, idx, ispin, bands)
cache_entries = [];
if ctx.sig.isdf.reuse_eps_real_wfn && ...
        isfield(ctx.eps, 'isdf_real_wfn') && ...
        size(ctx.eps.isdf_real_wfn, 1) >= wfn.ikq
    cache_entries = reshape(ctx.eps.isdf_real_wfn(wfn.ikq, ispin, :), ...
        1, ctx.nspinor);
end
components = isdf.real_components(wfn, fft_template, idx, ispin, ...
    ctx.nspinor, bands, cache_entries, false);
end

function tf = local_uses_vn_exchange(ctx)
exchange_space = lower(char(ctx.sig.isdf.exchange_space));
if ~any(strcmp(exchange_space, {'nn', 'vn'}))
    error('ISDF:UnknownExchangeSpace', ...
        'sig.isdf.exchange_space must be ''nn'' or ''vn'' (got ''%s'').', ...
        exchange_space);
end
tf = strcmp(exchange_space, 'vn') && ...
    ~ctx.sig.isdf.reuse_nn_for_vn;
end

function work = local_progress_work(block)
if isfield(block, 'progress') && isfield(block.progress, 'block_work')
    work = block.progress.block_work;
else
    work = 1;
end
end
