function polar = isdf_reduced_polarizability(vc_space, ev_occ, ev_unocc, options)
%ISDF_REDUCED_POLARIZABILITY Build reduced static polarizability coefficient.

if nargin < 4 || isempty(options)
    options = struct();
end

if isfield(vc_space, 'left_mu_components') && ...
        isfield(vc_space, 'right_mu_components')
    left_mu = vc_space.left_mu_components;
    right_mu = vc_space.right_mu_components;
elseif isfield(vc_space, 'phi_mu') && isfield(vc_space, 'psi_mu')
    % Scalar product spaces store the already-conjugated left factor.
    left_mu = {conj(vc_space.phi_mu)};
    right_mu = {vc_space.psi_mu};
else
    error('ISDF:ReducedPolarizabilitySpace', ...
        'vc_space must contain selected scalar or component wavefunctions.');
end

[coeff, info] = isdf_comega_cstar(left_mu, right_mu, ...
    ev_occ, ev_unocc, options);

polar = struct();
polar.coeff = coeff;
polar.info = info;
end
