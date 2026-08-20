function polar = polarizability(space, ev_occ, ev_unocc, options)
%ISDF.POLARIZABILITY Build reduced polarizability coefficients.

if nargin < 4 || isempty(options)
    options = struct();
end
if isfield(space, 'left_mu_components') && ...
        isfield(space, 'right_mu_components')
    left_mu = space.left_mu_components;
    right_mu = space.right_mu_components;
elseif isfield(space, 'phi_mu') && isfield(space, 'psi_mu')
    left_mu = {conj(space.phi_mu)};
    right_mu = {space.psi_mu};
else
    error('ISDF:ReducedPolarizabilitySpace', ...
        'space must contain selected scalar or component wavefunctions.');
end

[coeff, info] = comega_cstar( ...
    left_mu, right_mu, ev_occ, ev_unocc, options);
polar = struct('coeff', coeff, 'info', info);
end
