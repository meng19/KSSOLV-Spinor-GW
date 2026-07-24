script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

[sys, options, syms] = read_qe_gw( ...
    '.\example\qe_data\mos2_222_spinor', 1);
[sys, options] = gw_setup(sys, options);

eps.cutoff = 2;
eps.freq_dep = 0;
eps.nfreq = 1;
eps.freq = 0;
eps.inv = cell(sys.nkpts, 1);

sig.nbnd = 29;
sig.ndiag_min = 29;
sig.ndiag_max = 29;
sig.freq_dep = 0;
sig.coul_cut = 'spherical_truncation';
sig.coul_cutoff = 2;
sig.use_gpu = 0;
sig.precompute_wav = 0;
sig = sigma_set_defaults(sig);

probe = sigma_context(eps, sig, sys, options, syms, true);
for iq = 1:sys.nkpts
    eps.inv{iq} = eye(probe.sig.nmtx(iq));
end
ctx = sigma_context(eps, sig, sys, options, syms, false);

assert(strcmp(ctx.method, 'direct'));
assert(ctx.nk == sys.nkpts);
assert(ctx.nspinor == 2);
assert(numel(ctx.kdata) == sys.nkpts);
block = sigma_prepare_block(ctx, 1, 1, 29, 1, []);
assert(block.ik == 1 && block.iq == 1 && block.in == 29);
assert(isfield(block, 'coulg') && isfield(block, 'coulg_cutoff'));
assert(isfield(block, 'eps_inv'));
assert(numel(block.occ_kq) == sig.nbnd);

sig_reduced = sig;
sig_reduced.isdf.enable = true;
sig_reduced.isdf.algorithm = 'reduced_basis';
eps_empty = eps;
eps_empty.inv = cell(sys.nkpts, 1);
context_error_id = '';
try
    ctx_empty = sigma_context( ...
        eps_empty, sig_reduced, sys, options, syms, false);
catch ME
    context_error_id = ME.identifier;
end
assert(isempty(context_error_id), ...
    'Reduced context rejected a present inverse container with "%s".', ...
    context_error_id);
assert(all(cellfun(@isempty, ctx_empty.eps_inv_fbz)));

caught_id = '';
try
    sigma(eps_empty, sig_reduced, sys, options, syms);
catch ME
    caught_id = ME.identifier;
end
assert(strcmp(caught_id, 'ISDF:ReducedSigmaMissingQPoint'), ...
    'Expected per-q missing-screening error, got "%s".', caught_id);

star_fixture = tempname;
mkdir(star_fixture);
copyfile(fullfile(script_dir, 'fixtures', 'sigma_star_mismatch', ...
    'rqstar.txt'), fullfile(star_fixture, 'rqstar.m'));
original_path = path;
fixture_cleanup = onCleanup(@() rmdir(star_fixture, 's'));
path_cleanup = onCleanup(@() path(original_path));
addpath(star_fixture, '-begin');
clear rqstar;
caught_star_id = '';
try
    sigma(eps, sig_reduced, sys, options, syms);
catch ME
    caught_star_id = ME.identifier;
end
assert(strcmp(caught_star_id, 'ISDF:ReducedSigmaStar'), ...
    'Expected reduced q-star mismatch error, got "%s".', caught_star_id);
clear path_cleanup;
clear rqstar;
clear fixture_cleanup;

fprintf('ISDF sigma context/block test passed.\n');
