function polar = isdf_cauchy_polarizability(vc_space, ev_occ, ev_unocc, options)
%ISDF_CAUCHY_POLARIZABILITY Build reduced static polarizability coefficient.

if nargin < 4 || isempty(options)
    options = struct();
end

[coeff, info] = isdf_comega_cstar(vc_space.phi_mu, vc_space.psi_mu, ...
    ev_occ, ev_unocc, options);

polar = struct();
polar.coeff = coeff;
polar.info = info;
end
