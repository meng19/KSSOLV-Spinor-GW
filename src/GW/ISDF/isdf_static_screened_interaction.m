function screened = isdf_static_screened_interaction(vc_space, vcoul, polar)
%ISDF_STATIC_SCREENED_INTERACTION Construct reduced static screened operator.

vcoul = vcoul(:);
if size(vc_space.zeta_g, 1) ~= numel(vcoul)
    error('ISDF:ScreenedInteractionSize', ...
        'Coulomb vector length must match zeta_g row count.');
end

vmat = vc_space.zeta_g' * (vcoul .* vc_space.zeta_g);

warning_state_near = warning('query', 'MATLAB:nearlySingularMatrix');
warning_state_sing = warning('query', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
cleanup_warning = onCleanup(@() local_restore_warning(warning_state_near, warning_state_sing));

eps_mu = inv(polar.coeff) - vmat;
eps_mu_inv = (eye(size(polar.coeff)) - polar.coeff * vmat) \ polar.coeff;

screened = struct();
screened.zeta_g = vc_space.zeta_g;
screened.vcoul = vcoul;
screened.coeff = polar.coeff;
screened.vmat = vmat;
screened.eps_mu = eps_mu;
screened.eps_mu_inv = eps_mu_inv;
screened.w_mu = eps_mu_inv;
end

function local_restore_warning(warning_state_near, warning_state_sing)
warning(warning_state_near.state, 'MATLAB:nearlySingularMatrix');
warning(warning_state_sing.state, 'MATLAB:singularMatrix');
end
