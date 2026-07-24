function screened = screened_w(space, epsilon_vcoul, polar)
%ISDF.SCREENED_W Construct a reduced static screened operator.

epsilon_vcoul = epsilon_vcoul(:);
if size(space.zeta_g, 1) ~= numel(epsilon_vcoul)
    error('ISDF:ScreenedInteractionSize', ...
        'Coulomb vector length must match zeta_g row count.');
end
vmat = space.zeta_g' * (epsilon_vcoul .* space.zeta_g);

warning_state_near = warning('query', 'MATLAB:nearlySingularMatrix');
warning_state_sing = warning('query', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
cleanup_warning = onCleanup(@() local_restore_warning( ...
    warning_state_near, warning_state_sing)); %#ok<NASGU>

smw_denominator = inv(polar.coeff) - vmat;
k_mu = (eye(size(polar.coeff)) - polar.coeff * vmat) \ polar.coeff;
screened = struct();
screened.zeta_g = space.zeta_g;
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
