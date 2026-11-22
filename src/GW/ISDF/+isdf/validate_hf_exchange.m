function report = validate_hf_exchange(space, left, right, idx_q, fftgrid, coulg, pair_weights, options)
%VALIDATE_HF_EXCHANGE Compare direct and ISDF bare-exchange energies.
%
% The returned energies are positive Coulomb contractions, i.e. the
% magnitude before the conventional minus sign and GW prefactor are
% applied.  LEFT/RIGHT may contain one entry per spinor component.  RIGHT
% may instead be a provider struct with fields fun (a function handle that
% returns right components for a requested band list), bands (right band
% indices), and block_size (number of bands generated per streamed block).

if nargin < 8 || isempty(options)
    options = struct();
end
if ~iscell(left)
    left = {left};
end
provider = [];
if isstruct(right) && isfield(right, 'fun')
    % Streamed right components: regenerate right bands in bounded blocks
    % instead of receiving the full component set.
    provider = right;
    if ~isa(provider.fun, 'function_handle')
        error('ISDF:HFValidationProvider', ...
            'right.fun must be a function handle returning right components.');
    end
    right_bands = provider.bands(:).';
    nright = numel(right_bands);
    if ~(isnumeric(provider.block_size) && isscalar(provider.block_size) && ...
            provider.block_size >= 1)
        error('ISDF:HFValidationBlockSize', ...
            'right.block_size must be a positive scalar.');
    end
    block_size = min(nright, max(1, ceil(provider.block_size)));
else
    if ~iscell(right)
        right = {right};
    end
    if numel(left) ~= numel(right)
        error('ISDF:HFValidationComponents', ...
            'left and right must have the same number of spinor components.');
    end
    nright = size(right{1}, 2);
end
options = set_defaults(options, size(left{1}, 2), nright, ...
    size(left{1}, 1), fftgrid);
nleft = size(left{1}, 2);
npairs = nleft * nright;
if numel(pair_weights) ~= npairs
    error('ISDF:HFValidationWeights', ...
        'pair_weights must contain one value per left-right pair.');
end
if ~isfield(space, 'nleft') || ~isfield(space, 'nright') || ...
        space.nleft * space.nright ~= npairs
    error('ISDF:HFValidationSpace', ...
        'ISDF coefficients do not match the supplied product space.');
end
if size(space.zeta_g, 1) ~= numel(idx_q) || numel(coulg) ~= numel(idx_q)
    error('ISDF:HFValidationCoulomb', ...
        'idx_q, coulg, and space.zeta_g must have the same G dimension.');
end

max_pairs = options.validate_hf_exchange_max_pairs;
if ~(isnumeric(max_pairs) && isscalar(max_pairs) && max_pairs > 0)
    error('ISDF:HFValidationMaxPairs', ...
        'validate_hf_exchange_max_pairs must be a positive scalar.');
end
selected = find(pair_weights(:) ~= 0);
if isfinite(max_pairs) && numel(selected) > max_pairs
    selected = selected(1:floor(max_pairs));
end
weights = pair_weights(selected).';
label = options.validate_hf_exchange_label;

coulg = coulg(:)/1103.71209453158;
idx_q = idx_q(:);
direct_values = zeros(1, numel(selected));
isdf_values = zeros(1, numel(selected));
normalization = prod(fftgrid); % Must match private/zeta_to_g.m.
% Match the reference ISDF HF validator explicitly: c_rho' * tildeVq *
% c_rho.  This is algebraically identical to sum_G v_G |Z*c_rho|^2,
% but retains the reduced Coulomb matrix needed for like-for-like timing
% and numerical comparison with the original repository.
tilde_vq = space.zeta_g' * (coulg .* space.zeta_g);
if isempty(provider)
    for ipair = 1:numel(selected)
        column = selected(ipair);
        direct_values(ipair) = local_direct_pair_value( ...
            left, right, column, nleft, 1, idx_q, fftgrid, ...
            normalization, coulg);
        c_rho = isdf.product_coefficients(space, column);
        isdf_values(ipair) = gather(real(c_rho' * tilde_vq * c_rho));
    end
else
    % Streamed right components: evaluate every selected pair while its
    % right band is resident in the current block.
    for first = 1:block_size:nright
        last = min(nright, first + block_size - 1);
        right_block = provider.fun(right_bands(first:last));
        if ~iscell(right_block)
            right_block = {right_block};
        end
        if numel(right_block) ~= numel(left)
            error('ISDF:HFValidationComponents', ...
                ['Streamed right components must match the left ' ...
                'spinor count.']);
        end
        positions = find(selected >= (first - 1) * nleft + 1 & ...
            selected <= last * nleft);
        for jj = 1:numel(positions)
            ipair = positions(jj);
            column = selected(ipair);
            direct_values(ipair) = local_direct_pair_value( ...
                left, right_block, column, nleft, first, idx_q, ...
                fftgrid, normalization, coulg);
            c_rho = isdf.product_coefficients(space, column);
            isdf_values(ipair) = gather(real(c_rho' * tilde_vq * c_rho));
        end
    end
end

direct_energy = sum(weights .* direct_values);
isdf_energy = sum(weights .* isdf_values);
report = struct('direct_energy', direct_energy, ...
    'isdf_energy', isdf_energy, 'absolute_error', isdf_energy - direct_energy, ...
    'relative_error', abs(isdf_energy - direct_energy) / ...
        max(abs(direct_energy), eps), ...
    'pairs_total', nnz(pair_weights), 'pairs_checked', numel(selected), ...
    'max_pair_absolute_error', max(abs(isdf_values - direct_values)), ...
    'direct_pair_values', direct_values, 'isdf_pair_values', isdf_values, ...
    'pair_indices', selected, 'tilde_vq', tilde_vq);
fprintf(['ISDF %s validation: direct = %.12e, ISDF = %.12e, ' ...
    'abs err = %.3e, rel err = %.3e (%d/%d weighted pairs)\n'], ...
    label, report.direct_energy, report.isdf_energy, abs(report.absolute_error), ...
    report.relative_error, report.pairs_checked, report.pairs_total);
end

function value = local_direct_pair_value(left, right, column, nleft, ...
        first_band, idx_q, fftgrid, normalization, coulg)
% Direct Coulomb value of one left-right product pair.  RIGHT holds the
% spinor components of a contiguous band block whose first band is
% FIRST_BAND; COLUMN indexes the pair with the left band fastest.
ileft = mod(column - 1, nleft) + 1;
local_column = floor((column - 1) / nleft) + 1 - first_band + 1;
product = zeros(size(left{1}, 1), 1, 'like', left{1});
for icomponent = 1:numel(left)
    product = product + conj(left{icomponent}(:, ileft)) .* ...
        right{icomponent}(:, local_column);
end
product_g = fftn(reshape(product, fftgrid)) / normalization;
product_g = product_g(idx_q);
value = gather(real(sum(coulg .* abs(product_g).^2)));
end
