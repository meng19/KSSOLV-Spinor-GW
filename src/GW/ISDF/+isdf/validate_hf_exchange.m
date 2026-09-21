function report = validate_hf_exchange(space, left, right, idx_q, fftgrid, coulg, pair_weights, options)
%VALIDATE_HF_EXCHANGE Compare direct and ISDF bare-exchange energies.
%
% The returned energies are positive Coulomb contractions, i.e. the
% magnitude before the conventional minus sign and GW prefactor are
% applied.  LEFT/RIGHT may contain one entry per spinor component.

if nargin < 8 || isempty(options)
    options = struct();
end
if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
if numel(left) ~= numel(right)
    error('ISDF:HFValidationComponents', ...
        'left and right must have the same number of spinor components.');
end
nleft = size(left{1}, 2);
nright = size(right{1}, 2);
npairs = nleft * nright;
if numel(pair_weights) ~= npairs
    error('ISDF:HFValidationWeights', ...
        'pair_weights must contain one value per left-right pair.');
end
if size(space.product_mu, 2) ~= npairs
    error('ISDF:HFValidationSpace', ...
        'ISDF coefficients do not match the supplied product space.');
end
if size(space.zeta_g, 1) ~= numel(idx_q) || numel(coulg) ~= numel(idx_q)
    error('ISDF:HFValidationCoulomb', ...
        'idx_q, coulg, and space.zeta_g must have the same G dimension.');
end

max_pairs = local_option(options, 'validate_hf_exchange_max_pairs', Inf);
if ~(isnumeric(max_pairs) && isscalar(max_pairs) && max_pairs > 0)
    error('ISDF:HFValidationMaxPairs', ...
        'validate_hf_exchange_max_pairs must be a positive scalar.');
end
selected = find(pair_weights(:) ~= 0);
if isfinite(max_pairs) && numel(selected) > max_pairs
    selected = selected(1:floor(max_pairs));
end
weights = pair_weights(selected).';
label = local_option(options, 'validate_hf_exchange_label', 'HF exchange');

coulg = coulg(:)/1103.71209453158;
idx_q = idx_q(:);
direct_values = zeros(1, numel(selected));
isdf_values = zeros(1, numel(selected));
normalization = 1103.71209453158; % Must match private/zeta_to_g.m.
% Match the reference ISDF HF validator explicitly: c_rho' * tildeVq *
% c_rho.  This is algebraically identical to sum_G v_G |Z*c_rho|^2,
% but retains the reduced Coulomb matrix needed for like-for-like timing
% and numerical comparison with the original repository.
tilde_vq = space.zeta_g' * (coulg .* space.zeta_g);
for ipair = 1:numel(selected)
    column = selected(ipair);
    ileft = mod(column - 1, nleft) + 1;
    iright = floor((column - 1) / nleft) + 1;
    product = zeros(size(left{1}, 1), 1, 'like', left{1});
    for icomponent = 1:numel(left)
        product = product + conj(left{icomponent}(:, ileft)) .* ...
            right{icomponent}(:, iright);
    end
    product_g = ifftn(reshape(product, fftgrid)) * normalization;
    product_g = product_g(idx_q);
    direct_values(ipair) = gather(real(sum(coulg .* abs(product_g).^2)));
    c_rho = space.product_mu(:, column);
    isdf_values(ipair) = gather(real(c_rho' * tilde_vq * c_rho));
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

function value = local_option(options, name, default_value)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = default_value;
end
end
