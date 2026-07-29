function screened = screened_w(space, epsilon_vcoul, polar)
%ISDF.SCREENED_W Construct a reduced screened operator.

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

nmu = size(polar.coeff, 1);
nfreq = size(polar.coeff, 3);
k_mu = zeros(nmu, nmu, nfreq);
identity = eye(nmu);
for ifreq = 1:nfreq
    coeff = polar.coeff(:, :, ifreq);
    system_matrix = identity - coeff * vmat;
    k_mu(:, :, ifreq) = system_matrix \ coeff;
end
screened = struct();
screened.zeta_g = space.zeta_g;
screened.epsilon_vcoul = epsilon_vcoul;
screened.polar_coeff = polar.coeff;
screened.coulomb_reduced = vmat;
screened.k_mu = k_mu;
end

% ---- Warning-state cleanup ----

function local_restore_warning(warning_state_near, warning_state_sing)
warning(warning_state_near.state, 'MATLAB:nearlySingularMatrix');
warning(warning_state_sing.state, 'MATLAB:singularMatrix');
end
