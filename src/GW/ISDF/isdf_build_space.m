function space = isdf_build_space(phi, psi, idx_q, fftgrid, options)
%ISDF_BUILD_SPACE Build a compact ISDF product-space representation.

if nargin < 5 || isempty(options)
    options = struct();
end

options = isdf_set_defaults(options, size(phi, 2), size(psi, 2), size(phi, 1));
if ~isfield(options, 'fftgrid') || isempty(options.fftgrid)
    options.fftgrid = fftgrid;
end
ind_mu = isdf_indices(phi, psi, options);
[zeta_g, solve_info] = isdf_kernelg_current_fft(phi, psi, ind_mu, idx_q, fftgrid, options);

space = struct();
space.phi = phi;
space.psi = psi;
space.ind_mu = ind_mu;
space.zeta_g = zeta_g;
space.phi_mu = phi(ind_mu, :);
space.psi_mu = psi(ind_mu, :);
space.rank = numel(ind_mu);
space.options = options;
space.solve_info = solve_info;
end
