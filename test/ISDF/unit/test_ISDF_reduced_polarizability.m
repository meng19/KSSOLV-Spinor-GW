script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(12);
nmu = 5;
nv = 3;
nc = 4;
vc_space = struct();
vc_space.phi_mu = randn(nmu, nv) + 1i * randn(nmu, nv);
vc_space.psi_mu = randn(nmu, nc) + 1i * randn(nmu, nc);

ev_occ = [-0.8; -0.4; -0.1];
ev_unocc = [0.3; 0.7; 1.0; 1.4];

direct_options = struct('method', 'direct');
direct = isdf_reduced_polarizability(vc_space, ev_occ, ev_unocc, direct_options);

cauchy_options = struct('method', 'cauchy', 'froErr', 1e-10, 'MaxIter', 8);
cauchy = isdf_reduced_polarizability(vc_space, ev_occ, ev_unocc, cauchy_options);

relerr = norm(direct.coeff - cauchy.coeff, 'fro') / max(1, norm(direct.coeff, 'fro'));
assert(relerr < 1e-8, 'Reduced polarizability differs from direct: %.3e', relerr);

component_space = struct();
component_space.left_mu_components = { ...
    randn(nmu, nv) + 1i * randn(nmu, nv), ...
    randn(nmu, nv) + 1i * randn(nmu, nv)};
component_space.right_mu_components = { ...
    randn(nmu, nc) + 1i * randn(nmu, nc), ...
    randn(nmu, nc) + 1i * randn(nmu, nc)};
component_polar = isdf_reduced_polarizability( ...
    component_space, ev_occ, ev_unocc, direct_options);
component_ref = isdf_comega_cstar( ...
    component_space.left_mu_components, ...
    component_space.right_mu_components, ...
    ev_occ, ev_unocc, direct_options);
assert(norm(component_polar.coeff - component_ref, 'fro') < 1e-12, ...
    'Reduced polarizability must preserve spinor component cross terms.');

fprintf('ISDF reduced polarizability test passed. relerr = %.3e\n', relerr);
