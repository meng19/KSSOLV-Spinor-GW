function kernel = screened_kernel(screened, target_zeta_g, contract_vcoul)
%ISDF.SCREENED_KERNEL Reduced kernel for (epsilon^{-1}-I)*v.

if isfield(screened, 'epsilon_vcoul')
    epsilon_vcoul = screened.epsilon_vcoul;
elseif isfield(screened, 'vcoul')
    epsilon_vcoul = screened.vcoul;
else
    error('ISDF:ScreenedKernelMissingCoulomb', ...
        'screened must contain epsilon_vcoul.');
end
if nargin < 3 || isempty(contract_vcoul)
    contract_vcoul = epsilon_vcoul;
end
epsilon_vcoul = epsilon_vcoul(:);
contract_vcoul = contract_vcoul(:);
build_full_matrix = isempty(target_zeta_g);
if (~build_full_matrix && ...
        size(target_zeta_g, 1) ~= numel(epsilon_vcoul)) || ...
        size(screened.zeta_g, 1) ~= numel(epsilon_vcoul) || ...
        numel(contract_vcoul) ~= numel(epsilon_vcoul)
    error('ISDF:ScreenedKernelSize', ...
        ['target_zeta_g and Coulomb vectors must have matching ' ...
         'G dimensions.']);
end

if build_full_matrix
    left_projector = epsilon_vcoul .* screened.zeta_g;
    if isequal(epsilon_vcoul, contract_vcoul) && isreal(epsilon_vcoul)
        right_projector = left_projector';
    else
        right_projector = screened.zeta_g' .* contract_vcoul.';
    end
else
    left_projector = target_zeta_g.' * ...
        (epsilon_vcoul .* screened.zeta_g);
    if isequal(epsilon_vcoul, contract_vcoul) && isreal(epsilon_vcoul)
        right_projector = left_projector';
    else
        right_projector = screened.zeta_g' * ...
            (contract_vcoul .* conj(target_zeta_g));
    end
end
if ndims(screened.k_mu) == 2
    kernel = left_projector * screened.k_mu * right_projector;
else
    kernel = zeros(size(left_projector, 1), size(right_projector, 2), ...
        size(screened.k_mu, 3));
    for ifreq = 1:size(screened.k_mu, 3)
        kernel(:, :, ifreq) = left_projector * ...
            screened.k_mu(:, :, ifreq) * right_projector;
    end
end
end
