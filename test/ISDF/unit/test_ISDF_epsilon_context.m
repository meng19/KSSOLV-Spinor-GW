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

eps.nbnd = 29;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
eps.freq_dep = 0;
eps.cutoff = 2;
eps.coul_cut = 'spherical_truncation';
eps.coul_cutoff = 2;
eps.use_gpu = 0;
eps.precompute_wav = 0;
eps = epsilon_set_defaults(eps);

ctx = epsilon_context(sys, options, syms, eps);
assert(strcmp(ctx.method, 'direct'));
assert(ctx.nq == sys.nkpts);
assert(ctx.nspinor == 2);
assert(numel(ctx.qdata) == sys.nkpts);
assert(ctx.qdata{1}.nrq >= 1);
assert(numel(ctx.qdata{1}.g_maps) == ctx.qdata{1}.nrq);

block = epsilon_prepare_block(ctx, 1, 1, 1, []);
assert(block.iq == 1 && block.ik == 1 && block.ispin == 1);
assert(isfield(block, 'wfnk') && isfield(block, 'wfnkq'));
assert(isfield(block, 'fft') && isfield(block, 'idx'));
assert(numel(block.valence_bands) == sum(block.occ_vkq > 0));
assert(all(block.conduction_bands > sum(block.occ_ck > 0)));

fprintf('ISDF epsilon context/block test passed.\n');
