function screened = isdf_static_screened_interaction(vc_space, epsilon_vcoul, polar)
%ISDF_STATIC_SCREENED_INTERACTION Construct reduced static screened operator.

epsilon_vcoul = epsilon_vcoul(:);
if size(vc_space.zeta_g, 1) ~= numel(epsilon_vcoul)
    error('ISDF:ScreenedInteractionSize', ...
        'Coulomb vector length must match zeta_g row count.');
end

% epsilon_vcoul is the bare Coulomb used in epsilon = I - v * chi0.
vmat = vc_space.zeta_g' * (epsilon_vcoul .* vc_space.zeta_g);

warning_state_near = warning('query', 'MATLAB:nearlySingularMatrix');
warning_state_sing = warning('query', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
cleanup_warning = onCleanup(@() local_restore_warning(warning_state_near, warning_state_sing));

smw_denominator = inv(polar.coeff) - vmat;
k_mu = (eye(size(polar.coeff)) - polar.coeff * vmat) \ polar.coeff;

screened = struct();
screened.zeta_g = vc_space.zeta_g;
screened.epsilon_vcoul = epsilon_vcoul;
screened.polar_coeff = polar.coeff;
screened.coulomb_reduced = vmat;
screened.smw_denominator = smw_denominator;
screened.k_mu = k_mu;
end

function local_restore_warning(warning_state_near, warning_state_sing)
warning(warning_state_near.state, 'MATLAB:nearlySingularMatrix');
warning(warning_state_sing.state, 'MATLAB:singularMatrix');
end
