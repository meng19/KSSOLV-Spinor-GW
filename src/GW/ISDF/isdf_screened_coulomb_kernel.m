function kernel = isdf_screened_coulomb_kernel(screened, target_zeta_g, right_vcoul)
%ISDF_SCREENED_COULOMB_KERNEL Reduced kernel for (epsilon^{-1}-I)*v.

if nargin < 3 || isempty(right_vcoul)
    right_vcoul = screened.vcoul;
end

right_vcoul = right_vcoul(:);
left_vcoul = screened.vcoul(:);
if size(target_zeta_g, 1) ~= numel(left_vcoul) || numel(right_vcoul) ~= numel(left_vcoul)
    error('ISDF:ScreenedKernelSize', ...
        'target_zeta_g and Coulomb vectors must have matching G dimensions.');
end

left_projector = target_zeta_g.' * (left_vcoul .* screened.zeta_g);
right_projector = screened.zeta_g' * (right_vcoul .* conj(target_zeta_g));
kernel = left_projector * screened.eps_mu_inv * right_projector;
end
