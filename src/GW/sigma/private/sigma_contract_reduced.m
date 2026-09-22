function contribution = sigma_contract_reduced(ctx, block, matrix_elements)
%SIGMA_CONTRACT_REDUCED Reduced-basis screened-interaction contraction.

if isempty(block.screened_w)
    if isempty(block.eps_inv)
        error('ISDF:ReducedSigmaMissingQPoint', ...
            'No screened interaction is available for full-BZ q-point %d.', ...
            block.iq_fbz);
    end
    contribution = sigma_contract_full(ctx, block, matrix_elements);
    return;
end

if size(block.screened_w.zeta_g, 1) ~= block.n_cutoff
    error('ISDF:ReducedSigmaScreenedSize', ...
        'Reduced screened interaction does not match sigma cutoff.');
end

gw_timer('start', 'Sigma screened kernel');
kernel = local_project_screened_kernel(ctx, block, matrix_elements.space, ...
    local_space_key(matrix_elements));
gw_timer('stop', 'Sigma screened kernel');
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
if ctx.sig.freq_dep == 2
    omega = [];
    iw_lda = [];
    % Bare exchange is frequency independent. Compute all occupied-band
    % terms once in the same ISDF coefficient-space metric used below.
    occ = reshape(block.occ_kq, 1, []);
    ax_loc = local_batch_exchange(ctx, block, matrix_elements, occ);
end
progress_work = gw_block_work(block, 1);
if ctx.sig.freq_dep == 0
    gw_timer('start', 'Sigma static contraction');
    [asx_loc, ax_loc, ach_loc] = local_static_batch_contract( ...
        ctx, block, matrix_elements, kernel(:, :, 1));
    gw_timer('stop', 'Sigma static contraction');
    gw_block_progress(block, progress_work, ...
        sprintf('S b%d k%d q%d static: all %d sum bands', ...
        block.in, block.ik, block.iq, ctx.nbands));
else
if ctx.sig.freq_dep ~= 2
    error('ISDF:ReducedSigmaFrequencyMode', ...
        'Reduced sigma supports freq_dep = 0 or 2 (got %g).', ...
        ctx.sig.freq_dep);
end
for nn = 1:ctx.nbands
    if isfield(matrix_elements, 'coeff')
        coeff = matrix_elements.coeff(:, nn);
    else
        coeff = matrix_elements.space.product_mu(:, nn);
    end
    [asx_loc, ach_loc, achx_loc_nn(block.in, nn), ...
        omega, iw_lda] = sigma_fullfreq(asx_loc, ach_loc, ...
        block.in, nn, block.wfnk.ikq, block.wfnkq.ikq, ...
        block.occ_kq(nn), ctx.options.ev, block.ispin, ...
        coeff, coeff, ctx.fact * kernel, ctx.sig);
    gw_block_progress(block, progress_work * (0.5 + 0.5 * nn / ctx.nbands), ...
        sprintf('S b%d i%d q%d n%d/%d', ...
        block.in, block.ik, block.iq, nn, ctx.nbands));
end
end

achx_loc = 0;
if ctx.sig.exact_static_ch
    if ctx.sig.freq_dep == 0
        % Use the established exact-CH path: form the G-space screened
        % interaction and then apply the precomputed G-G' difference map.
        % This dense route trades an nG-by-nG temporary for faster BLAS.
        gw_timer('start', 'Sigma full screened kernel');
        screened_matrix = ctx.fact * local_full_kernel(ctx, block);
        gw_timer('stop', 'Sigma full screened kernel');
        exact_ch = sigma_cohsex_exact_ch(block.in, block.ispin, ...
            ctx.fbz, ctx.kdata{block.ik}.indrk, block.iq, block.aqsch, ...
            screened_matrix, ctx.sig, block.igpp, block.valid_indices);
        achx_loc = sum(exact_ch, 'all');
    elseif ctx.sig.freq_dep == 2
        gw_timer('start', 'Sigma full screened kernel');
        screened_matrix = ctx.fact * local_full_kernel(ctx, block);
        gw_timer('stop', 'Sigma full screened kernel');
        kdata = ctx.kdata{block.ik};
        exact_ch = sigma_cohsex_exact_ch(block.in, block.ispin, ...
            ctx.fbz, kdata.indrk, block.iq, block.aqsch, ...
            screened_matrix, ctx.sig, block.igpp, block.valid_indices);
        achx_loc_nn(block.in, 1) = achx_loc_nn(block.in, 1) + ...
            0.5 * 0.5 * sum(exact_ch, 'all');
        achx_loc = sum(achx_loc_nn(block.in, :), 'all');
    end
end

if ctx.sig.freq_dep == 2
    contribution = sigma_make_contribution( ...
        ctx, asx_loc, ax_loc, ach_loc, achx_loc, ...
        omega, iw_lda, achx_loc_nn);
else
    contribution = sigma_make_contribution( ...
        ctx, asx_loc, ax_loc, ach_loc, achx_loc);
end
end

function [asx_loc, ax_loc, ach_loc] = local_static_batch_contract( ...
        ctx, block, matrix_elements, kernel)
% Evaluate all static NN contractions as dense matrix products.  For C
% containing one ISDF coefficient column per summation band, the desired
% values are diag(C.' * kernel * conj(C)).

coeff = local_coefficients(matrix_elements, ctx.nbands);
kernel_coeff = kernel * conj(coeff);
screened_values = ctx.fact * sum(coeff .* kernel_coeff, 1);

occ = reshape(block.occ_kq, 1, []);
if isa(screened_values, 'gpuArray') && ~isa(occ, 'gpuArray')
    occ = gpuArray(occ);
end
if local_has_vn_screened_exchange(matrix_elements)
    vn = matrix_elements.gme_exchange;
    gw_timer('start', 'Sigma VN screened kernel');
    vn_kernel = local_project_screened_kernel(ctx, block, vn.space, ...
        vn.space_key);
    gw_timer('stop', 'Sigma VN screened kernel');
    vn_kernel_coeff = vn_kernel * conj(vn.coeff);
    vn_values = ctx.fact * sum(vn.coeff .* vn_kernel_coeff, 1);
    asx_loc = -sum(occ(vn.bands) .* vn_values);
else
    asx_loc = -sum(occ .* screened_values);
end
ax_loc = local_batch_exchange(ctx, block, matrix_elements, occ);
ach_loc = sum(screened_values);
end

function tf = local_has_vn_screened_exchange(matrix_elements)
tf = isfield(matrix_elements, 'gme_exchange') && ...
    isstruct(matrix_elements.gme_exchange) && ...
    isfield(matrix_elements.gme_exchange, 'space') && ...
    isfield(matrix_elements.gme_exchange, 'coeff') && ...
    ~isempty(matrix_elements.gme_exchange.space) && ...
    ~isempty(matrix_elements.gme_exchange.coeff);
end

function ax_loc = local_batch_exchange(ctx, block, matrix_elements, occ)
if isfield(matrix_elements, 'gme_exchange') && ...
        isstruct(matrix_elements.gme_exchange)
    vn = matrix_elements.gme_exchange;
    bands = vn.bands;
    if isfield(vn, 'space') && isfield(vn, 'coeff') && ...
            ~isempty(vn.space) && ~isempty(vn.coeff)
        gw_timer('start', 'Sigma VN bare exchange GME');
        exchange_values = local_bare_exchange_values( ...
            block, vn.space, vn.coeff);
        gw_timer('stop', 'Sigma VN bare exchange GME');
        ax_loc = -ctx.fact * sum(occ(bands) .* exchange_values);
        return;
    elseif isfield(vn, 'values') && ~isempty(vn.values)
        values = vn.values;
        exchange_values = sum(bsxfun(@times, abs(values).^2, ...
            block.coulg), 1);
        ax_loc = -ctx.fact * sum(occ(bands) .* exchange_values);
        return;
    end
end
if isfield(matrix_elements, 'space') && isfield(matrix_elements, 'coeff') && ...
        ~isempty(matrix_elements.space) && ~isempty(matrix_elements.coeff)
    % Bare exchange has support only on occupied bands.  Do not form GME
    % columns for empty states that will be multiplied by zero afterwards.
    occupied = find(gw_gather_if_gpu(occ) > 0);
    if isempty(occupied)
        ax_loc = 0;
        return;
    end
    gw_timer('start', 'Sigma NN bare exchange GME');
    exchange_values = local_bare_exchange_values( ...
        block, matrix_elements.space, matrix_elements.coeff(:, occupied));
    gw_timer('stop', 'Sigma NN bare exchange GME');
    ax_loc = -ctx.fact * sum(occ(occupied) .* exchange_values);
    return;
end
exchange_values = sum(bsxfun(@times, abs(matrix_elements.gme).^2, ...
    block.coulg), 1);
ax_loc = -ctx.fact * sum(occ .* exchange_values);
end

function coeff = local_coefficients(matrix_elements, nbands)
if isfield(matrix_elements, 'coeff')
    coeff = matrix_elements.coeff;
else
    coeff = matrix_elements.space.product_mu(:, 1:nbands);
end
end

function kernel = local_project_screened_kernel(ctx, block, space, space_key)
% Project W-v into one target ISDF product space.  It is shared by NN and
% VN contractions; their distinct space keys keep cached projections safe.
target_zeta = space.zeta_g(1:block.n_cutoff, :);
kernel = isdf.screened_kernel(block.screened_w, target_zeta, ...
    block.coulg_cutoff, local_kernel_key(ctx, block, 'target', ...
    space_key));
end

function values = local_bare_exchange_values(block, space, coeff)
% G-space bare exchange: sum_G v_G*abs((zeta_g*coeff)_Gn)^2.
target_zeta = space.zeta_g(1:numel(block.coulg), :);
gme = target_zeta * coeff;
values = sum(bsxfun(@times, abs(gme).^2, block.coulg), 1);
end

function kernel = local_full_kernel(ctx, block)
% q-only projection used by the exact static Coulomb-hole reference.

kernel = isdf.screened_kernel(block.screened_w, [], block.coulg_cutoff, ...
    local_kernel_key(ctx, block, 'full', ''));
end

function key = local_space_key(matrix_elements)
key = '';
if isstruct(matrix_elements) && isfield(matrix_elements, 'space_key') && ...
        ischar(matrix_elements.space_key)
    key = matrix_elements.space_key;
end
end

function key = local_kernel_key(ctx, block, kind, space_key)
% Reuse is safe only when the caller can name the target product space and
% the configured byte budget is positive; otherwise return '' so that the
% kernel is recomputed per block.

key = '';
if strcmp(kind, 'target') && isempty(space_key)
    return;
end
if isfield(ctx.sig.isdf, 'screened_kernel_cache_gb')
    budget = ctx.sig.isdf.screened_kernel_cache_gb;
elseif isfield(ctx.sig.isdf, 'screened_kernel_cache_bytes')
    budget = ctx.sig.isdf.screened_kernel_cache_bytes / (1024^3);
else
    return;
end
if ~(isnumeric(budget) && isscalar(budget) && budget > 0)
    return;
end

% The key must distinguish every input of ISDF.SCREENED_KERNEL: the target
% space (space_key), the full-BZ q-point that fixes screened_w and the
% cutoff-limited Coulomb vector (iq_fbz, n_cutoff), the reduced dimensions
% and the CPU/GPU kind of the stored arrays.
zeta_g = block.screened_w.zeta_g;
k_mu = block.screened_w.k_mu;
key = sprintf('%s|%s|q%d|n%d|v%d|m%d|p%d|g%d', kind, space_key, ...
    block.iq_fbz, block.n_cutoff, size(zeta_g, 2), size(k_mu, 1), ...
    size(k_mu, 3), double(isa(zeta_g, 'gpuArray')));
end
