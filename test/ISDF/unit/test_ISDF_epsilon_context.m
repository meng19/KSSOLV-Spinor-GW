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

eps_reduced = eps;
eps_reduced.isdf.enable = true;
eps_reduced.isdf.algorithm = 'reduced_basis';
ctx_reduced = epsilon_context(sys, options, syms, eps_reduced);

ctx_invalid = ctx_reduced;
ctx_invalid.eps.isdf.output = 'invalid_output';
invalid_output_id = local_ops_error_id(ctx_invalid);
assert(strcmp(invalid_output_id, 'ISDF:ReducedEpsilonOutput'), ...
    'Expected output validation during ops construction, got "%s".', ...
    invalid_output_id);

acc_full = local_reduced_accumulator(ctx_reduced, 'full_inverse');
assert(acc_full.need_full_inverse && ~acc_full.need_screened_w);
assert(isequal(size(acc_full.chi0), ...
    [ctx_reduced.pol.nmtx(1), ctx_reduced.pol.nmtx(1)]));
acc_screened = local_reduced_accumulator(ctx_reduced, 'screened_w');
assert(~acc_screened.need_full_inverse && acc_screened.need_screened_w);
assert(isempty(acc_screened.chi0));
acc_both = local_reduced_accumulator(ctx_reduced, 'both');
assert(acc_both.need_full_inverse && acc_both.need_screened_w);
assert(isequal(size(acc_both.chi0), ...
    [ctx_reduced.pol.nmtx(1), ctx_reduced.pol.nmtx(1)]));

star_fixture = tempname;
mkdir(star_fixture);
copyfile(fullfile(script_dir, 'fixtures', 'epsilon_star_mismatch', ...
    'rqstar.txt'), fullfile(star_fixture, 'rqstar.m'));
original_path = path;
fixture_cleanup = onCleanup(@() local_cleanup_fixture( ...
    original_path, star_fixture));
addpath(star_fixture, '-begin');
clear rqstar;

direct_star_id = local_epsilon_error_id( ...
    sys, options, syms, eps);
assert(strcmp(direct_star_id, 'MATLAB:error:nonScalarInput'), ...
    'Expected legacy direct k-star error, got "%s".', direct_star_id);

eps_matrix = eps;
eps_matrix.isdf.enable = true;
eps_matrix.isdf.algorithm = 'matrix_elements';
matrix_star_id = local_epsilon_error_id( ...
    sys, options, syms, eps_matrix);
assert(strcmp(matrix_star_id, 'MATLAB:error:nonScalarInput'), ...
    'Expected legacy matrix-elements k-star error, got "%s".', ...
    matrix_star_id);

reduced_star_id = local_epsilon_error_id( ...
    sys, options, syms, eps_reduced);
assert(strcmp(reduced_star_id, 'ISDF:ReducedEpsilonStar'), ...
    'Expected reduced k-star error, got "%s".', reduced_star_id);

clear fixture_cleanup;

fprintf('ISDF epsilon context/block test passed.\n');

% ---- Error-capture helpers ----

function identifier = local_epsilon_error_id(sys, options, syms, eps)
identifier = '';
try
    epsilon(sys, options, syms, eps);
catch ME
    identifier = ME.identifier;
end
end

function identifier = local_ops_error_id(ctx)
identifier = '';
try
    epsilon_ops(ctx);
catch ME
    identifier = ME.identifier;
end
end

% ---- Reduced accumulator construction ----

function acc = local_reduced_accumulator(ctx, output)
ctx.eps.isdf.output = output;
ops = epsilon_ops(ctx);
acc = ops.init(1);
end

% ---- Temporary path cleanup ----

function local_cleanup_fixture(original_path, star_fixture)
path(original_path);
clear rqstar;
rmdir(star_fixture, 's');
end
