# ISDF Shared-Loop Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor epsilon and sigma so direct, ISDF matrix-elements, and ISDF reduced-basis calculations share their physical main loops while algorithm-specific numerical operations are selected once through context/block/ops interfaces.

**Architecture:** `epsilon.m` retains one q/spin/k loop and `sigma.m` retains one spin/band/k/q loop. Immutable context structs hold reusable setup, block structs hold one loop position, strategy function tables provide algorithm-specific evaluation/contraction, and ISDF numerical kernels live under the MATLAB `+isdf` package.

**Tech Stack:** MATLAB, KSSOLV GW data structures, MATLAB package folders, script-based MATLAB tests, Git.

## Global Constraints

- Start execution from an isolated worktree based on the branch containing this plan; the worktree must include commit `2bc7108`.
- The recoverable pre-refactor implementation is tag `isdf-before-shared-loop-refactor` at commit `e186aab`.
- Keep public signatures `epsilon(sys, options, syms, eps)` and `sigma(eps, sig, sys, options, syms)` unchanged.
- Keep existing `eps.isdf.*`, `sig.isdf.*`, `eps.inv`, `eps.isdf_screened_w`, `sig.sig`, `sig.cor`, and `sig.eqp0` semantics.
- Reduced-basis remains static and CPU-only; full-frequency and GPU support are not expanded.
- `eps.isdf.output='screened_w'` must not construct or store full `eps.inv`.
- Main loops may call preselected `ops` handles but must not parse `isdf.enable` or `isdf.algorithm` inside the loops.
- Preserve `eps.save_mem`, `eps.precompute_wav`, `sig.precompute_wav`, `sig.no_symmetries_q_grid`, and `sig.exact_static_ch` behavior.
- Do not weaken existing numerical tolerances.
- Use `apply_patch` for source changes and run `git diff --check` before every commit.

---

## Locked File Structure

### Workflow files

```text
src/GW/common/gw_resolve_method.m

src/GW/epsilon/epsilon.m
src/GW/epsilon/epsilon_context.m
src/GW/epsilon/epsilon_prepare_block.m
src/GW/epsilon/epsilon_ops.m

src/GW/sigma/sigma.m
src/GW/sigma/sigma_context.m
src/GW/sigma/sigma_prepare_block.m
src/GW/sigma/sigma_ops.m
```

### Package files

```text
src/GW/ISDF/+isdf/build_space.m
src/GW/ISDF/+isdf/polarizability.m
src/GW/ISDF/+isdf/screened_w.m
src/GW/ISDF/+isdf/screened_kernel.m
src/GW/ISDF/+isdf/screened_contract.m
src/GW/ISDF/+isdf/real_component.m
src/GW/ISDF/+isdf/matrix_elements.m

src/GW/ISDF/+isdf/private/product_gram.m
src/GW/ISDF/+isdf/private/component_products.m
src/GW/ISDF/+isdf/private/sample_points.m
src/GW/ISDF/+isdf/private/stable_solve.m
src/GW/ISDF/+isdf/private/zeta_to_g.m
```

### Locked interfaces

```matlab
method = gw_resolve_method(isdf_options, workflow)

ctx = epsilon_context(sys, options, syms, eps)
block = epsilon_prepare_block(ctx, iq, ik, ispin, prepared)
ops = epsilon_ops(ctx)

ctx = sigma_context(eps, sig, sys, options, syms, metadata_only)
block = sigma_prepare_block(ctx, ik, iq, in, ispin, prepared)
ops = sigma_ops(ctx)

space = isdf.build_space(left_components, right_components, idx_q, fftgrid, options)
polar = isdf.polarizability(space, ev_occ, ev_unocc, options)
screened = isdf.screened_w(space, epsilon_vcoul, polar)
kernel = isdf.screened_kernel(screened, target_zeta_g, contract_vcoul)
value = isdf.screened_contract(kernel, coeff_left, coeff_right)
values = isdf.real_component(wfn, fft_template, idx, ispin, ispinor, bands)
gme = isdf.matrix_elements(left_components, right_components, idx_q, fftgrid, options)
```

`left_components` always stores unconjugated physical wavefunctions. Every package product uses

```matlab
P(r,ij) = sum_s conj(left_components{s}(r,i)) * right_components{s}(r,j)
```

This convention removes the legacy ambiguity where scalar callers sometimes passed an already conjugated `phi`.

---

### Task 1: Canonical Method Resolution

**Files:**
- Create: `src/GW/common/gw_resolve_method.m`
- Create: `test/ISDF/unit/test_ISDF_method_resolution.m`

**Interfaces:**
- Consumes: an ISDF settings struct and workflow name `'epsilon'` or `'sigma'`.
- Produces: canonical char vector `'direct'`, `'matrix_elements'`, or `'reduced_basis'`.

- [ ] **Step 1: Write the failing method-resolution test**

Create `test/ISDF/unit/test_ISDF_method_resolution.m`:

```matlab
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

assert(strcmp(gw_resolve_method(struct('enable', false), 'epsilon'), 'direct'));
assert(strcmp(gw_resolve_method(struct( ...
    'enable', true, 'algorithm', 'matrix_elements'), 'epsilon'), ...
    'matrix_elements'));
assert(strcmp(gw_resolve_method(struct( ...
    'enable', true, 'algorithm', 'REDUCED_BASIS'), 'sigma'), ...
    'reduced_basis'));

try
    gw_resolve_method(struct('enable', true, 'algorithm', 'bad'), 'epsilon');
    error('TEST:ExpectedFailure', 'Unknown epsilon method did not fail.');
catch ME
    assert(strcmp(ME.identifier, 'ISDF:UnknownEpsilonAlgorithm'));
end

try
    gw_resolve_method(struct('enable', true, 'algorithm', 'bad'), 'sigma');
    error('TEST:ExpectedFailure', 'Unknown sigma method did not fail.');
catch ME
    assert(strcmp(ME.identifier, 'ISDF:UnknownSigmaAlgorithm'));
end

fprintf('ISDF method resolution test passed.\n');
```

- [ ] **Step 2: Run the test and verify the missing-function failure**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_method_resolution.m')"
```

Expected: FAIL with `Undefined function 'gw_resolve_method'`.

- [ ] **Step 3: Implement the canonical resolver**

Create `src/GW/common/gw_resolve_method.m`:

```matlab
function method = gw_resolve_method(isdf_options, workflow)
%GW_RESOLVE_METHOD Resolve one canonical GW workflow method.

workflow = lower(char(workflow));
if ~any(strcmp(workflow, {'epsilon', 'sigma'}))
    error('ISDF:UnknownWorkflow', 'Unknown GW workflow "%s".', workflow);
end

if nargin < 1 || isempty(isdf_options) || ~isstruct(isdf_options) || ...
        ~isfield(isdf_options, 'enable') || ~isdf_options.enable
    method = 'direct';
    return;
end
if ~isfield(isdf_options, 'algorithm') || isempty(isdf_options.algorithm)
    error_id = sprintf('ISDF:Unknown%sAlgorithm', ...
        [upper(workflow(1)), workflow(2:end)]);
    error(error_id, '%s ISDF algorithm is missing.', workflow);
end

method = lower(char(isdf_options.algorithm));
if ~any(strcmp(method, {'matrix_elements', 'reduced_basis'}))
    error_id = sprintf('ISDF:Unknown%sAlgorithm', ...
        [upper(workflow(1)), workflow(2:end)]);
    error(error_id, ...
        ['Unknown %s ISDF algorithm "%s". Supported algorithms: ' ...
        'reduced_basis, matrix_elements.'], workflow, method);
end
end
```

- [ ] **Step 4: Run the test and existing defaults guard**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_method_resolution.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_cauchy_cohsex_guards.m')"
```

Expected: both PASS.

- [ ] **Step 5: Commit the method contract**

```powershell
git add -- src/GW/common/gw_resolve_method.m test/ISDF/unit/test_ISDF_method_resolution.m
git diff --cached --check
git commit -m "Refactor GW ISDF method resolution"
```

---

### Task 2: Introduce the `+isdf` Public Numerical API

**Files:**
- Create: `src/GW/ISDF/+isdf/build_space.m`
- Create: `src/GW/ISDF/+isdf/polarizability.m`
- Create: `src/GW/ISDF/+isdf/screened_w.m`
- Create: `src/GW/ISDF/+isdf/screened_kernel.m`
- Create: `src/GW/ISDF/+isdf/screened_contract.m`
- Create: `src/GW/ISDF/+isdf/real_component.m`
- Create: `src/GW/ISDF/+isdf/matrix_elements.m`
- Create: `test/ISDF/unit/test_ISDF_package_api.m`

**Interfaces:**
- Consumes: physical left/right component cells using the locked unconjugated-left convention.
- Produces: namespaced package calls usable by epsilon/sigma strategies while legacy implementations remain the numerical reference.

- [ ] **Step 1: Write the failing package API test**

Create `test/ISDF/unit/test_ISDF_package_api.m`:

```matlab
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

rng(61, 'twister');
fftgrid = [4, 3, 2];
ngrid = prod(fftgrid);
idx_q = [1; 2; 5; 9; 13; 24];
left = {randn(ngrid, 2) + 1i * randn(ngrid, 2), ...
    randn(ngrid, 2) + 1i * randn(ngrid, 2)};
right = {randn(ngrid, 3) + 1i * randn(ngrid, 3), ...
    randn(ngrid, 3) + 1i * randn(ngrid, 3)};
options = struct('rank', 6, 'sample_method', 'qrcp', 'seed', 0);

space = isdf.build_space(left, right, idx_q, fftgrid, options);
gme = isdf.matrix_elements(left, right, idx_q, fftgrid, options);
reconstructed = space.zeta_g * space.product_mu;
assert(max(abs(gme(:) - reconstructed(:))) < 1e-12);

legacy_space = isdf_build_space(left, right, idx_q, fftgrid, options);
assert(max(abs(space.zeta_g(:) - legacy_space.zeta_g(:))) < 1e-12);
assert(isequal(space.ind_mu, legacy_space.ind_mu));

ev_occ = [-0.8; -0.2];
ev_unocc = [0.3; 0.9; 1.4];
solver = struct('method', 'direct');
polar = isdf.polarizability(space, ev_occ, ev_unocc, solver);
legacy_polar = isdf_reduced_polarizability( ...
    legacy_space, ev_occ, ev_unocc, solver);
assert(norm(polar.coeff - legacy_polar.coeff, 'fro') < 1e-12);

vcoul = 0.5 + rand(numel(idx_q), 1);
screened = isdf.screened_w(space, vcoul, polar);
legacy_screened = isdf_static_screened_interaction(space, vcoul, polar);
assert(norm(screened.k_mu - legacy_screened.k_mu, 'fro') < 1e-12);

target = randn(numel(idx_q), 4) + 1i * randn(numel(idx_q), 4);
kernel = isdf.screened_kernel(screened, target, vcoul);
legacy_kernel = isdf_screened_coulomb_kernel( ...
    legacy_screened, target, vcoul);
assert(norm(kernel - legacy_kernel, 'fro') < 1e-12);

coeff = randn(4, 1) + 1i * randn(4, 1);
assert(abs(isdf.screened_contract(kernel, coeff) - ...
    isdf_screened_coulomb_contract(legacy_kernel, coeff)) < 1e-12);

fprintf('ISDF package API compatibility test passed.\n');
```

- [ ] **Step 2: Run the package test and verify namespace failure**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_package_api.m')"
```

Expected: FAIL because `isdf.build_space` is undefined.

- [ ] **Step 3: Add compatibility package functions**

Create the package functions with these exact bodies:

```matlab
% +isdf/build_space.m
function space = build_space(left, right, idx_q, fftgrid, options)
if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
space = isdf_build_space(left, right, idx_q, fftgrid, options);
end

% +isdf/polarizability.m
function polar = polarizability(space, ev_occ, ev_unocc, options)
polar = isdf_reduced_polarizability(space, ev_occ, ev_unocc, options);
end

% +isdf/screened_w.m
function screened = screened_w(space, epsilon_vcoul, polar)
screened = isdf_static_screened_interaction(space, epsilon_vcoul, polar);
end

% +isdf/screened_kernel.m
function kernel = screened_kernel(screened, target_zeta_g, contract_vcoul)
kernel = isdf_screened_coulomb_kernel( ...
    screened, target_zeta_g, contract_vcoul);
end

% +isdf/screened_contract.m
function value = screened_contract(kernel, coeff_left, coeff_right)
if nargin < 3
    coeff_right = [];
end
value = isdf_screened_coulomb_contract( ...
    kernel, coeff_left, coeff_right);
end

% +isdf/real_component.m
function values = real_component(wfn, fft_template, idx, ispin, ispinor, bands)
values = isdf_wavefunction_real_component( ...
    wfn, fft_template, idx, ispin, ispinor, bands);
end
```

Create `src/GW/ISDF/+isdf/matrix_elements.m`:

```matlab
function gme = matrix_elements(left, right, idx_q, fftgrid, options)
%ISDF.MATRIX_ELEMENTS Fourier matrix elements for component products.

if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
space = isdf.build_space(left, right, idx_q, fftgrid, options);
nleft = size(left{1}, 2);
nright = size(right{1}, 2);
gme = reshape(space.zeta_g * space.product_mu, ...
    numel(idx_q), nleft, nright);
end
```

- [ ] **Step 4: Run package and existing numerical unit tests**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_package_api.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_component_product_space.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_static_screened_interaction.m')"
```

Expected: all PASS.

- [ ] **Step 5: Commit the package façade**

```powershell
git add -- src/GW/ISDF/+isdf test/ISDF/unit/test_ISDF_package_api.m
git diff --cached --check
git commit -m "Add namespaced ISDF numerical API"
```

---

### Task 3: Extract Epsilon Context and Block Preparation

**Files:**
- Create: `src/GW/epsilon/epsilon_context.m`
- Create: `src/GW/epsilon/epsilon_prepare_block.m`
- Create: `test/ISDF/unit/test_ISDF_epsilon_context.m`
- Modify: `src/GW/epsilon/epsilon.m:1-127`

**Interfaces:**
- Consumes: defaulted epsilon input and optional prepared wavefunction/FFT data.
- Produces: immutable context plus one normalized q/k/spin block.

- [ ] **Step 1: Write the failing epsilon context test**

Create `test/ISDF/unit/test_ISDF_epsilon_context.m`:

```matlab
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
```

- [ ] **Step 2: Run the test and verify missing context failure**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_epsilon_context.m')"
```

Expected: FAIL with `Undefined function 'epsilon_context'`.

- [ ] **Step 3: Implement `epsilon_context`**

Move the common setup from `epsilon.m:2-77` and the star-map logic from `epsilon.m:143-165` into `epsilon_context.m`. The resulting function must expose these fields:

```matlab
function ctx = epsilon_context(sys, options, syms, eps)
ctx.sys = sys;
ctx.options = options;
ctx.syms = syms;
ctx.eps = eps;
ctx.method = gw_resolve_method(eps.isdf, 'epsilon');
ctx.ryd = 13.6056923;
ctx.nbands = eps.nbnd;
ctx.nspin = sys.nspin;
ctx.nspinor = sys.nspinor;
ctx.nq = sys.nkpts;
ctx.wfc_cutoff = 2 * sys.ecut;
ctx.save_mem = eps.save_mem;
ctx.precompute_wav = eps.precompute_wav;
ctx.use_gpu = eps.use_gpu && exist('gpuDevice', 'file');

if strcmp(ctx.method, 'matrix_elements') && ctx.use_gpu
    error('ISDF:EpsilonGPUUnsupported', ...
        'ISDF epsilon currently supports CPU execution only.');
end
if strcmp(ctx.method, 'reduced_basis') && ctx.use_gpu
    error('ISDF:ReducedEpsilonGPU', ...
        'ISDF reduced-basis epsilon currently supports CPU execution only.');
end
if strcmp(ctx.method, 'reduced_basis') && eps.freq_dep ~= 0
    error('ISDF:ReducedEpsilonFrequency', ...
        'ISDF reduced-basis epsilon requires eps.freq_dep = 0.');
end

sigrid = Ggrid(sys, 4 * sys.ecut);
ctx.gvec = Gvector(sigrid, sys);
ctx.gr = fullbz(options, syms, true);
ctx.fact = 4 / (ctx.gr.nf * sys.vol * ctx.nspin * ctx.nspinor);
ctx.pol.qpt = options.kpts;

if eps.freq_dep == 2 && eps.freq_dep_method == 2
    ctx.pol.nfreq_rel = fix(eps.freq_cutoff / eps.delta_freq) + 1;
    ctx.pol.nfreq = ctx.pol.nfreq_rel + eps.nfreq_imag;
    ctx.pol.freq_grid = 0:eps.delta_freq:eps.freq_cutoff;
    tmp = 0:(eps.nfreq_imag - 1);
    ctx.pol.freq_brd = -2i * ctx.ryd * tmp ./ (tmp - eps.nfreq_imag);
    ctx.pol.freq = [ctx.pol.freq_grid, ctx.pol.freq_brd];
else
    ctx.pol.nfreq = 1;
    ctx.pol.freq = 0;
end

ctx.ekin = zeros(ctx.gvec.ng, ctx.nq);
ctx.qdata = cell(ctx.nq, 1);
for iq = 1:ctx.nq
    qq = ctx.pol.qpt(iq, :);
    [ctx.ekin(:, iq), ctx.pol.isrtx(:, iq)] = sortrx( ...
        qq, ctx.gvec.ng, ctx.gvec.mill, sys);
    ctx.pol.nmtx(:, iq) = gcutoff(ctx.gvec.ng, ctx.ekin(:, iq), ...
        ctx.pol.isrtx(:, iq), eps.cutoff);
    ctx.pol.mtx{:, iq} = ctx.gvec.mill( ...
        ctx.pol.isrtx(1:ctx.pol.nmtx(iq), iq), :);
    box_min = zeros(1, 3);
    box_max = zeros(1, 3);
    [box_min, box_max] = get_gvecs_bounds( ...
        ctx.pol.mtx{:, iq}, box_min, box_max);
    ctx.pol.fftgrid{:, iq} = min( ...
        options.wfn_fftgrid + box_max - box_min, options.fftgrid);
    ctx.qdata{iq} = local_qdata(ctx, iq);
end
end
```

Add this local `local_qdata(ctx, iq)`; it preserves the existing star ordering while making the identity map explicit for star member 1:

```matlab
function qdata = local_qdata(ctx, iq)
qq = ctx.pol.qpt(iq, :);
syms_q = subgrp(qq, ctx.syms);
[nrq, neq, indrk] = irrbz(syms_q, ctx.gr);

isorti = zeros(ctx.gvec.ng, 1);
for ig = 1:ctx.gvec.ng
    isorti(ctx.pol.isrtx(ig, iq)) = ig;
end

qdata.nrq = nrq;
qdata.neq = neq;
qdata.indrk = indrk;
qdata.syms = syms_q;
qdata.g_maps = cell(nrq, 1);
qdata.rqs = cell(nrq, 1);
for ik = 1:nrq
    rk = ctx.gr.f(indrk(ik), :);
    [nstar, indst, rqs] = rqstar(syms_q, rk);
    if nstar ~= neq(ik)
        error('ISDF:ReducedEpsilonStar', ...
            'K-point star size does not match its irreducible weight.');
    end
    qdata.rqs{ik} = rqs;
    qdata.g_maps{ik} = cell(nstar, 1);
    qdata.g_maps{ik}{1} = (1:ctx.pol.nmtx(iq)).';
    for it = 2:nstar
        itran = syms_q.indsub(indst(it));
        kgq = -syms_q.kgzero(indst(it), :);
        qdata.g_maps{ik}{it} = gmap( ...
            ctx.gvec, ctx.syms, ctx.pol.nmtx(iq), itran, kgq, ...
            ctx.pol.isrtx(:, iq), isorti, ctx.sys);
    end
end
end
```

- [ ] **Step 4: Implement `epsilon_prepare_block`**

Create `src/GW/epsilon/epsilon_prepare_block.m`:

```matlab
function block = epsilon_prepare_block(ctx, iq, ik, ispin, prepared)
qdata = ctx.qdata{iq};
qq = ctx.pol.qpt(iq, :);
rk = ctx.gr.f(qdata.indrk(ik), :);

if nargin >= 5 && ~isempty(prepared)
    wfnk = prepared.wfnk;
    wfnkq = prepared.wfnkq;
    fft = prepared.fft;
    idx = prepared.idx;
else
    wfnk = genwf(rk, ctx.gr, ctx.gvec, ctx.syms, ctx.sys, ...
        ctx.options, ctx.wfc_cutoff, ctx.nbands, ctx.use_gpu);
    wfnkq = genwf(rk + qq, ctx.gr, ctx.gvec, ctx.syms, ctx.sys, ...
        ctx.options, ctx.wfc_cutoff, ctx.nbands, ctx.use_gpu);
    [fft, idx] = epsilon_prefft( ...
        wfnkq, wfnk, iq, ik, ctx.pol, [], [], ctx.use_gpu);
end

block.iq = iq;
block.ik = ik;
block.ispin = ispin;
block.ik_fbz = qdata.indrk(ik);
block.q = qq;
block.k = rk;
block.weight = qdata.neq(ik);
block.g_maps = qdata.g_maps{ik};
block.star_kpoints = qdata.rqs{ik};
block.wfnk = wfnk;
block.wfnkq = wfnkq;
block.fft = fft;
block.idx = idx;
block.occ_vkq = get_occ(ctx.options, wfnkq.ikq, ispin);
block.occ_ck = get_occ(ctx.options, wfnk.ikq, ispin);
block.valence_bands = 1:sum(block.occ_vkq > 0);
block.conduction_bands = (sum(block.occ_ck > 0) + 1):ctx.nbands;
block.ev_occ = ctx.options.ev( ...
    block.valence_bands, wfnkq.ikq, ispin);
block.ev_unocc = ctx.options.ev( ...
    block.conduction_bands, wfnk.ikq, ispin);
end
```

- [ ] **Step 5: Replace duplicated setup reads in `epsilon.m` without changing its loops**

At the start of `epsilon.m`, after defaults, create `ctx` and bind the current local names from it:

```matlab
ctx = epsilon_context(sys, options, syms, eps);
ryd = ctx.ryd;
nbands = ctx.nbands;
nspin = ctx.nspin;
nspinor = ctx.nspinor;
wfc_cutoff = ctx.wfc_cutoff;
save_mem = ctx.save_mem;
use_gpu = ctx.use_gpu;
gvec = ctx.gvec;
gr = ctx.gr;
pol = ctx.pol;
ekin = ctx.ekin;
fact = ctx.fact;
```

Do not remove either legacy main path in this task. This step only makes setup output come from one context.

- [ ] **Step 6: Run context and epsilon regression tests**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_epsilon_context.m')"
matlab -wait -batch "run('test/ISDF/validation/test_mos2_222_isdf_epsilon_validation.m')"
```

Expected: both PASS; the matrix-elements validation keeps its existing tolerance.

- [ ] **Step 7: Commit epsilon context extraction**

```powershell
git add -- src/GW/epsilon/epsilon.m src/GW/epsilon/epsilon_context.m src/GW/epsilon/epsilon_prepare_block.m test/ISDF/unit/test_ISDF_epsilon_context.m
git diff --cached --check
git commit -m "Extract shared epsilon context and blocks"
```

---

### Task 4: Share the Epsilon Loop for Direct and Matrix-Elements Modes

**Files:**
- Create: `src/GW/epsilon/epsilon_ops.m`
- Create: `test/ISDF/unit/test_ISDF_epsilon_shared_loop.m`
- Modify: `src/GW/epsilon/epsilon.m:78-356`

**Interfaces:**
- Consumes: `ctx` and normalized epsilon blocks.
- Produces: one preselected ops table with `init`, `evaluate`, `accumulate`, and `finalize` handles.

- [ ] **Step 1: Write the failing structural test**

Create `test/ISDF/unit/test_ISDF_epsilon_shared_loop.m`:

```matlab
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
epsilon_file = fullfile(repo_root, 'src', 'GW', 'epsilon', 'epsilon.m');
source = fileread(epsilon_file);

assert(contains(source, 'ops = epsilon_ops(ctx);'));
assert(contains(source, 'contribution = ops.evaluate(block);'));
assert(contains(source, 'acc = ops.accumulate(acc, contribution, block);'));
assert(~contains(source, 'use_isdf = isfield'), ...
    'epsilon.m must not resolve ISDF mode inside its main loop.');

fprintf('ISDF shared epsilon loop structure test passed.\n');
```

- [ ] **Step 2: Run the structural test and verify failure**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_epsilon_shared_loop.m')"
```

Expected: FAIL because `epsilon.m` does not contain the ops calls.

- [ ] **Step 3: Implement direct and matrix-elements ops**

Create `src/GW/epsilon/epsilon_ops.m` with this public selector:

```matlab
function ops = epsilon_ops(ctx)
ops.name = ctx.method;
switch ctx.method
    case 'direct'
        ops.init = @(iq) local_init_full(ctx, iq);
        ops.evaluate = @(block) local_evaluate_full(ctx, block, false);
        ops.accumulate = @(acc, contribution, block) ...
            local_accumulate_full(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) local_finalize_full(ctx, eps, acc, iq);
    case 'matrix_elements'
        ops.init = @(iq) local_init_full(ctx, iq);
        ops.evaluate = @(block) local_evaluate_full(ctx, block, true);
        ops.accumulate = @(acc, contribution, block) ...
            local_accumulate_full(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) local_finalize_full(ctx, eps, acc, iq);
    case 'reduced_basis'
        error('ISDF:ReducedEpsilonNotIntegrated', ...
            'Reduced epsilon is integrated in the next migration task.');
end
end
```

`local_evaluate_full` must return one struct with `gme` and `eden`:

```matlab
contribution.gme = zeros(numel(block.idx.q), ...
    numel(block.valence_bands), numel(block.conduction_bands));
if use_isdf
    left = cell(1, ctx.nspinor);
    right = cell(1, ctx.nspinor);
    for ispinor = 1:ctx.nspinor
        left{ispinor} = isdf.real_component(block.wfnkq, ...
            block.fft.Nfft1, block.idx.kq, block.ispin, ispinor, ...
            block.valence_bands);
        right{ispinor} = isdf.real_component(block.wfnk, ...
            block.fft.Nfft2, block.idx.k, block.ispin, ispinor, ...
            block.conduction_bands);
    end
    isdf_options = ctx.eps.isdf;
    if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
        isdf_options.rank = ceil(sqrt(numel(block.valence_bands) * ...
            numel(block.conduction_bands)) * ctx.eps.isdf.rank_ratio);
    end
    contribution.gme = isdf.matrix_elements(left, right, ...
        block.idx.q, size(block.fft.Nfft1), isdf_options);
else
    for iv_local = 1:numel(block.valence_bands)
        iv = block.valence_bands(iv_local);
        for ic_local = 1:numel(block.conduction_bands)
            ic = block.conduction_bands(ic_local);
            contribution.gme(:, iv_local, ic_local) = getm_epsilon( ...
                iv, ic, block.wfnkq, block.wfnk, block.fft, block.idx, ...
                block.ispin, ctx.nspinor, ctx.use_gpu);
        end
    end
end
```

Fill `contribution.eden(iv_local, ic_local, ifreq)` using the current `get_eden` call from `epsilon.m:257-260` for every band pair and frequency.

`local_accumulate_full` must preserve both memory modes:

- when `ctx.save_mem` is true, immediately add each star-mapped matrix element outer product to `acc.chi0`;
- when false, append `{gme, eden}` entries to `acc.deferred{ispin, block.ik_fbz}` and defer the same outer products until finalize.

Use this exact helper for both paths:

```matlab
function chi0 = local_add_states(chi0, gme, eden)
for iv = 1:size(gme, 2)
    for ic = 1:size(gme, 3)
        vector = gme(:, iv, ic);
        for ifreq = 1:size(eden, 3)
            chi0(:, :, ifreq) = chi0(:, :, ifreq) + ...
                conj(vector) * vector.' * eden(iv, ic, ifreq);
        end
    end
end
end
```

`local_finalize_full` must evaluate deferred entries, multiply by `ctx.fact`, build Coulomb with `coulG_select`, invert each frequency page exactly as the current CPU/GPU code does, and store `eps.inv{iq}`.

- [ ] **Step 4: Rewrite the direct/matrix main loop around ops**

In `epsilon.m`, select ops once before `for iq`:

```matlab
ops = epsilon_ops(ctx);
```

For direct and matrix-elements modes, replace the duplicated per-band algorithm condition with:

```matlab
for iq = 1:ctx.nq
    qdata = ctx.qdata{iq};
    acc = ops.init(iq);
    for ispin = 1:ctx.nspin
        for ik = 1:qdata.nrq
            prepared = local_epsilon_prepared_data( ...
                ctx, iq, ik, wfnk_all, wfnkq_all, fft_all, idx_all);
            block = epsilon_prepare_block(ctx, iq, ik, ispin, prepared);
            if isempty(block.valence_bands) || isempty(block.conduction_bands)
                continue;
            end
            contribution = ops.evaluate(block);
            acc = ops.accumulate(acc, contribution, block);
        end
    end
    eps = ops.finalize(eps, acc, iq);
end
```

`local_epsilon_prepared_data` returns `[]` when precomputation is disabled; otherwise it returns the four fields expected by `epsilon_prepare_block`. Keep progress reporting outside strategy functions.

Leave the existing reduced-basis early dispatch temporarily in place; Task 5 removes it.

- [ ] **Step 5: Run direct and matrix-elements epsilon tests**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_epsilon_shared_loop.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_epsilon_matrix_elements.m')"
matlab -wait -batch "run('test/ISDF/validation/test_mos2_222_isdf_epsilon_validation.m')"
matlab -wait -batch "run('test/ISDF/validation/test_AgBr_isdf_epsilon_validation.m')"
```

Expected: all PASS at their existing tolerances.

- [ ] **Step 6: Commit the first shared epsilon loop**

```powershell
git add -- src/GW/epsilon/epsilon.m src/GW/epsilon/epsilon_ops.m test/ISDF/unit/test_ISDF_epsilon_shared_loop.m
git diff --cached --check
git commit -m "Share epsilon loop for direct and ISDF matrices"
```

---

### Task 5: Integrate Reduced-Basis Epsilon into the Shared Loop

**Files:**
- Modify: `src/GW/epsilon/epsilon.m:1-356`
- Modify: `src/GW/epsilon/epsilon_ops.m`
- Modify: `test/ISDF/unit/test_ISDF_epsilon_shared_loop.m`
- Delete: `src/GW/ISDF/isdf_epsilon_reduced_basis.m`

**Interfaces:**
- Consumes: the same epsilon block used by direct/matrix modes.
- Produces: low-rank contributions and final `screened_w`, `full_inverse`, or `both` output without a second workflow loop.

- [ ] **Step 1: Strengthen the structural test before removing the legacy path**

Append to `test_ISDF_epsilon_shared_loop.m`:

```matlab
assert(~contains(source, 'isdf_epsilon_reduced_basis'));
assert(isempty(regexp(source, ...
    'isdf_epsilon_reduced_basis\([^;]+;\s*return', 'once')));
assert(exist(fullfile(repo_root, 'src', 'GW', 'ISDF', ...
    'isdf_epsilon_reduced_basis.m'), 'file') == 0);
```

- [ ] **Step 2: Run the strengthened test and verify legacy-reference failure**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_epsilon_shared_loop.m')"
```

Expected: FAIL because the legacy file and early dispatch still exist.

- [ ] **Step 3: Add reduced epsilon ops**

Replace the reduced error branch in `epsilon_ops` with handles:

```matlab
case 'reduced_basis'
    ops.init = @(iq) local_init_reduced(ctx, iq);
    ops.evaluate = @(block) local_evaluate_reduced(ctx, block);
    ops.accumulate = @(acc, contribution, block) ...
        local_accumulate_reduced(ctx, acc, contribution, block);
    ops.finalize = @(eps, acc, iq) ...
        local_finalize_reduced(ctx, eps, acc, iq);
```

Implement `local_init_reduced` with the existing output contract, so screened-W-only mode never allocates a dense inverse accumulator:

```matlab
function acc = local_init_reduced(ctx, iq)
output_mode = lower(ctx.eps.isdf.output);
acc.need_full_inverse = any(strcmp(output_mode, {'full_inverse', 'both'}));
acc.need_screened_w = any(strcmp(output_mode, {'screened_w', 'both'}));
if ~acc.need_full_inverse && ~acc.need_screened_w
    error('ISDF:ReducedEpsilonOutput', ...
        'Unknown ISDF reduced-basis epsilon output "%s".', ...
        ctx.eps.isdf.output);
end
if acc.need_full_inverse
    acc.chi0 = zeros(ctx.pol.nmtx(iq));
else
    acc.chi0 = [];
end
acc.zeta_blocks = {};
acc.coeff_blocks = {};
acc.rank = cell(ctx.nspin, ctx.qdata{iq}.nrq);
acc.info = cell(ctx.nspin, ctx.qdata{iq}.nrq);
end
```

Implement `local_evaluate_reduced` by moving the component construction, rank default, `isdf.build_space`, and `isdf.polarizability` operations from `isdf_epsilon_reduced_basis.m:103-127`:

```matlab
left = cell(1, ctx.nspinor);
right = cell(1, ctx.nspinor);
for ispinor = 1:ctx.nspinor
    left{ispinor} = isdf.real_component(block.wfnkq, ...
        block.fft.Nfft1, block.idx.kq, block.ispin, ispinor, ...
        block.valence_bands);
    right{ispinor} = isdf.real_component(block.wfnk, ...
        block.fft.Nfft2, block.idx.k, block.ispin, ispinor, ...
        block.conduction_bands);
end
isdf_options = ctx.eps.isdf;
if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
    npairs = numel(block.valence_bands) * numel(block.conduction_bands);
    isdf_options.rank = ceil(sqrt(npairs) * ctx.eps.isdf.rank_ratio);
end
space = isdf.build_space(left, right, block.idx.q, ...
    size(block.fft.Nfft1), isdf_options);
solver.method = ctx.eps.isdf.reduced_solver;
solver.froErr = ctx.eps.isdf.cauchy_froErr;
solver.MaxIter = ctx.eps.isdf.cauchy_MaxIter;
polar = isdf.polarizability( ...
    space, block.ev_occ, block.ev_unocc, solver);
contribution.space = space;
contribution.polar = polar;
```

`local_accumulate_reduced` must apply every `block.g_maps{it}` exactly as `isdf_epsilon_reduced_basis.m:129-144`, appending `conj(zeta_star)` and `conj(polar.coeff)` to block lists. Store per-block metadata inside the accumulator without embedding the outer q index:

```matlab
acc.rank{block.ispin, block.ik} = contribution.space.rank;
acc.info{block.ispin, block.ik} = contribution.polar.info;
```

In `local_finalize_reduced`, copy that metadata to the public three-dimensional cells:

```matlab
for ispin = 1:size(acc.rank, 1)
    for ik = 1:size(acc.rank, 2)
        eps.isdf_reduced_rank{iq, ispin, ik} = acc.rank{ispin, ik};
        eps.isdf_reduced_info{iq, ispin, ik} = acc.info{ispin, ik};
    end
end
```

`local_finalize_reduced` must copy the output semantics from `isdf_epsilon_reduced_basis.m:146-159`:

```matlab
coulg = coulG_select(ctx.eps, ctx.pol.nmtx(iq), ...
    ctx.pol.isrtx(:, iq), ctx.ekin(:, iq), 0, ...
    ctx.pol.mtx{:, iq}, ctx.gvec, ctx.sys, iq);
if acc.need_full_inverse
    epsilon_matrix = eye(ctx.pol.nmtx(iq)) - ...
        coulg(:) .* (ctx.fact * acc.chi0);
    eps.inv{iq} = inv(epsilon_matrix);
end
if acc.need_screened_w && ~isempty(acc.zeta_blocks)
    combined_space.zeta_g = cat(2, acc.zeta_blocks{:});
    combined_polar.coeff = ctx.fact * blkdiag(acc.coeff_blocks{:});
    eps.isdf_screened_w{iq} = isdf.screened_w( ...
        combined_space, coulg(:), combined_polar);
end
```

- [ ] **Step 4: Remove the reduced early dispatch and legacy workflow file**

Delete the `use_isdf_reduced` block at `epsilon.m:22-28`. Let all methods enter the one ops-based loop. Remove `src/GW/ISDF/isdf_epsilon_reduced_basis.m` after all code has moved into context/ops.

At the end of `epsilon.m`, always assign common metadata from `ctx.pol`; only assign `eps.inv` through the strategy that produces it.

- [ ] **Step 5: Run reduced epsilon and structural tests**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_epsilon_shared_loop.m')"
matlab -wait -batch "run('test/ISDF/validation/test_AgBr_isdf_cauchy_epsilon_validation.m')"
matlab -wait -batch "run('test/ISDF/validation/test_mos2_222_isdf_reduced_validation.m')"
```

Expected: all PASS. The MoS2 epsilon relative error remains below `1e-8`; default screened-W mode has no `eps.inv`.

- [ ] **Step 6: Commit shared reduced epsilon**

```powershell
git add -- src/GW/epsilon/epsilon.m src/GW/epsilon/epsilon_ops.m test/ISDF/unit/test_ISDF_epsilon_shared_loop.m
git add -u -- src/GW/ISDF/isdf_epsilon_reduced_basis.m
git diff --cached --check
git commit -m "Integrate reduced epsilon into shared loop"
```

---

### Task 6: Extract Sigma Context and Block Preparation

**Files:**
- Create: `src/GW/sigma/sigma_context.m`
- Create: `src/GW/sigma/sigma_prepare_block.m`
- Create: `test/ISDF/unit/test_ISDF_sigma_context.m`
- Modify: `src/GW/sigma/sigma.m:1-179`

**Interfaces:**
- Consumes: epsilon screening results, defaulted sigma input, and optional prepared wavefunctions/indices.
- Produces: shared full-BZ screening maps and one normalized target-band/q block.

- [ ] **Step 1: Write the failing sigma context test**

Create `test/ISDF/unit/test_ISDF_sigma_context.m`:

```matlab
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

fprintf('ISDF sigma context/block test passed.\n');
```

The sixth `metadata_only` argument is explicitly part of the context constructor during testing. When true, context builds grids/mappings but skips dereferencing screening matrices.

- [ ] **Step 2: Run the test and verify missing context failure**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_sigma_context.m')"
```

Expected: FAIL with `Undefined function 'sigma_context'`.

- [ ] **Step 3: Implement `sigma_context`**

Create `sigma_context(eps, sig, sys, options, syms, metadata_only)` and default `metadata_only=false` when omitted. Move the common setup from `sigma.m:2-179` plus mapped screened-W logic from `isdf_sigma_reduced_basis.m:20-102` into it.

The context must expose:

```matlab
ctx.eps = eps;
ctx.sig = sig;
ctx.sys = sys;
ctx.options = options;
ctx.syms = syms;
ctx.method = gw_resolve_method(sig.isdf, 'sigma');
ctx.ryd = 13.6056923;
ctx.nbands = sig.nbnd;
ctx.band_range = sig.ndiag_min:sig.ndiag_max;
ctx.nspin = sys.nspin;
ctx.nspinor = sys.nspinor;
ctx.nk = sys.nkpts;
ctx.wfc_cutoff = 2 * sys.ecut;
ctx.use_gpu = sig.use_gpu && exist('gpuDevice', 'file');
ctx.gr = fullbz(options, syms, true);
ctx.fact = 1 / (ctx.gr.nf * sys.vol);
ctx.grid_size = [sys.n1, sys.n2, sys.n3];
```

Preserve the public-entry guard behavior before allocating grids:

```matlab
if ~strcmp(ctx.method, 'direct') && ctx.use_gpu
    error('ISDF:SigmaGPUUnsupported', ...
        'ISDF sigma path currently supports CPU execution only. Set sig.use_gpu = 0.');
end
if strcmp(ctx.method, 'reduced_basis') && ...
        (eps.freq_dep ~= 0 || sig.freq_dep ~= 0)
    error('ISDF:ReducedSigmaFrequency', ...
        'ISDF reduced-basis sigma requires static epsilon and sigma.');
end
```

Initialize the shared FFT templates in the context and preserve the current GPU-allocation fallback exactly:

```matlab
if ctx.use_gpu
    try
        ctx.fft.Nfft1 = gpuArray.zeros(ctx.grid_size);
        ctx.fft.Nfft2 = gpuArray.zeros(ctx.grid_size);
    catch
        warning('GPU memory insufficient for FFT grids. Falling back to CPU.');
        ctx.use_gpu = false;
        ctx.fft.Nfft1 = zeros(ctx.grid_size);
        ctx.fft.Nfft2 = zeros(ctx.grid_size);
    end
else
    ctx.fft.Nfft1 = zeros(ctx.grid_size);
    ctx.fft.Nfft2 = zeros(ctx.grid_size);
end
ctx.fft.size = prod(ctx.grid_size);
```

Validate method constraints before building loops. Build `ctx.sig.isrtx/nmtx/mtx`, `ctx.fbz`, `ctx.ekin_fbz`, `ctx.kdata`, `ctx.eps_inv_fbz`, and `ctx.screened_fbz`. Use one local `map_screened_w` with the exact row permutation and size check from `isdf_sigma_reduced_basis.m:259-270`.

For full-frequency direct/matrix modes, preserve all frequency fields currently initialized at `sigma.m:60-75`.

- [ ] **Step 4: Implement `sigma_prepare_block`**

Create `sigma_prepare_block(ctx, ik, iq, in, ispin, prepared)` so it:

1. reads `indrk/neq` from `ctx.kdata{ik}`;
2. obtains `wfnk/wfnkq/idx` from `prepared` or calls `genwf/sigma_prefft`;
3. computes occupations and Coulomb vectors;
4. attaches mapped `eps_inv` and/or `screened_w` for the full-BZ q index;
5. attaches exact-CH `igpp/valid_indices` from prepared data or `pre_exact_static_ch`.

The returned block fields are exactly:

```matlab
block.ik = ik;
block.iq = iq;
block.iq_fbz = iq_fbz;
block.in = in;
block.ispin = ispin;
block.weight = kdata.neq(iq);
block.wfnk = wfnk;
block.wfnkq = wfnkq;
block.idx = idx;
block.fft = ctx.fft;
block.occ_kq = get_occ(ctx.options, wfnkq.ikq, ispin);
block.n_cutoff = ctx.fbz.nmtx_cutoff(iq_fbz);
block.coulg = coulg;
block.coulg_cutoff = coulg(1:block.n_cutoff);
block.eps_inv = ctx.eps_inv_fbz{iq_fbz};
block.screened_w = ctx.screened_fbz{iq_fbz};
block.igpp = igpp;
block.valid_indices = valid_indices;
```

- [ ] **Step 5: Bind common sigma locals from context without changing the main loop**

After `sigma_set_defaults`, construct `ctx` and bind current locals from it. Keep the reduced early dispatch and both loops temporarily. The purpose of this task is one setup implementation, not loop migration.

- [ ] **Step 6: Run context and existing sigma unit tests**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_sigma_context.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_sigma_matrix_elements.m')"
```

Expected: both PASS.

- [ ] **Step 7: Commit sigma context extraction**

```powershell
git add -- src/GW/sigma/sigma.m src/GW/sigma/sigma_context.m src/GW/sigma/sigma_prepare_block.m test/ISDF/unit/test_ISDF_sigma_context.m
git diff --cached --check
git commit -m "Extract shared sigma context and blocks"
```

---

### Task 7: Share the Sigma Loop for Direct and Matrix-Elements Modes

**Files:**
- Create: `src/GW/sigma/sigma_ops.m`
- Create: `test/ISDF/unit/test_ISDF_sigma_shared_loop.m`
- Modify: `src/GW/sigma/sigma.m:180-379`

**Interfaces:**
- Consumes: one target-band/q block.
- Produces: matrix-element matrix plus a normalized contribution struct containing `asx`, `ax`, `ach`, `achx` and optional full-frequency fields.

- [ ] **Step 1: Write the failing structural test**

Create `test/ISDF/unit/test_ISDF_sigma_shared_loop.m`:

```matlab
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
sigma_file = fullfile(repo_root, 'src', 'GW', 'sigma', 'sigma.m');
source = fileread(sigma_file);

assert(contains(source, 'ops = sigma_ops(ctx);'));
assert(contains(source, 'matrix_elements = ops.matrix_elements(block);'));
assert(contains(source, 'contribution = ops.contract(block, matrix_elements);'));
assert(isempty(regexp(source, ...
    'if\s+use_isdf_matrix_elements', 'once')), ...
    'sigma.m must not select matrix construction inside the band loop.');

fprintf('ISDF shared sigma loop structure test passed.\n');
```

- [ ] **Step 2: Run the structural test and verify failure**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_sigma_shared_loop.m')"
```

Expected: FAIL because sigma has no ops calls.

- [ ] **Step 3: Implement direct and matrix-elements sigma ops**

Create `src/GW/sigma/sigma_ops.m`:

```matlab
function ops = sigma_ops(ctx)
ops.name = ctx.method;
switch ctx.method
    case 'direct'
        ops.matrix_elements = @(block) ...
            local_matrix_elements(ctx, block, false);
        ops.contract = @(block, matrix_elements) ...
            local_contract_full(ctx, block, matrix_elements);
    case 'matrix_elements'
        ops.matrix_elements = @(block) ...
            local_matrix_elements(ctx, block, true);
        ops.contract = @(block, matrix_elements) ...
            local_contract_full(ctx, block, matrix_elements);
    case 'reduced_basis'
        error('ISDF:ReducedSigmaNotIntegrated', ...
            'Reduced sigma is integrated in the next migration task.');
end
end
```

`local_matrix_elements` returns `[numel(block.idx.q), ctx.nbands]`. The direct branch loops `nn` and calls `getm_sigma`. The ISDF branch builds left/right component cells and calls the generic package API:

```matlab
left = cell(1, ctx.nspinor);
right = cell(1, ctx.nspinor);
for ispinor = 1:ctx.nspinor
    left{ispinor} = isdf.real_component(block.wfnk, ...
        block.fft.Nfft1, block.idx.k, block.ispin, ispinor, block.in);
    right{ispinor} = isdf.real_component(block.wfnkq, ...
        block.fft.Nfft2, block.idx.kq, block.ispin, ispinor, ...
        1:ctx.nbands);
end
isdf_options = ctx.sig.isdf;
if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
    isdf_options.rank = ceil(sqrt(ctx.nbands) * ctx.sig.isdf.rank_ratio);
end
gme3 = isdf.matrix_elements(left, right, block.idx.q, ...
    ctx.grid_size, isdf_options);
matrix_elements = reshape(gme3, numel(block.idx.q), ctx.nbands);
```

`local_contract_full` moves the existing `nn` loop from `sigma.m:293-321`. Return this exact common shape:

```matlab
contribution.asx = asx_loc;
contribution.ax = ax_loc;
contribution.ach = ach_loc;
contribution.achx = achx_loc;
contribution.omega = omega;
contribution.iw_lda = iw_lda;
contribution.asx_freq = asx_freq;
contribution.ach_freq = ach_freq;
contribution.achx_nn = achx_nn;
```

For static mode, unused frequency fields are empty arrays. Preserve all current `sigma_cohsex`, `sigma_fullfreq`, and exact-CH calls and factors.

- [ ] **Step 4: Rewrite the direct/matrix sigma loop around ops**

Select `ops = sigma_ops(ctx)` once before the spin loop. Replace only the algorithm-specific inner code with:

```matlab
block = sigma_prepare_block(ctx, ik, iq, in, ispin, prepared);
matrix_elements = ops.matrix_elements(block);
if ctx.sig.exact_static_ch && block.iq_fbz == 1
    aqsch{in, ispin} = matrix_elements(:, in);
end
block.aqsch = aqsch;
contribution = ops.contract(block, matrix_elements);
asxtemp = asxtemp + contribution.asx * block.weight;
axtemp = axtemp + contribution.ax * block.weight;
achtemp = achtemp + contribution.ach * block.weight;
if ctx.sig.exact_static_ch
    achxtemp = achxtemp + contribution.achx * block.weight;
end
```

`aqsch` is initialized once in `sigma.m` with the same lifetime as today. `local_contract_full` must read `block.aqsch` when it calls `sigma_cohsex_exact_ch`; do not hide this cross-q state in a persistent variable or closure.

Keep existing progress reporting, array storage, `quasi_energy`, and full-frequency finalization in `sigma.m`. Leave the reduced early dispatch until Task 8.

- [ ] **Step 5: Run sigma matrix and validation tests**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_sigma_shared_loop.m')"
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_sigma_matrix_elements.m')"
matlab -wait -batch "run('test/ISDF/validation/test_AgBr_isdf_sigma_validation.m')"
```

Expected: all PASS at existing tolerances.

- [ ] **Step 6: Commit the first shared sigma loop**

```powershell
git add -- src/GW/sigma/sigma.m src/GW/sigma/sigma_ops.m test/ISDF/unit/test_ISDF_sigma_shared_loop.m
git diff --cached --check
git commit -m "Share sigma loop for direct and ISDF matrices"
```

---

### Task 8: Integrate Reduced-Basis Sigma into the Shared Loop

**Files:**
- Modify: `src/GW/sigma/sigma.m:1-379`
- Modify: `src/GW/sigma/sigma_ops.m`
- Modify: `test/ISDF/unit/test_ISDF_sigma_shared_loop.m`
- Delete: `src/GW/ISDF/isdf_sigma_reduced_basis.m`

**Interfaces:**
- Consumes: the shared sigma block with either mapped screened-W or full inverse screening.
- Produces: the same contribution struct and final sigma arrays as direct/matrix modes.

- [ ] **Step 1: Strengthen the structural test**

Append:

```matlab
assert(~contains(source, 'isdf_sigma_reduced_basis'));
assert(isempty(regexp(source, ...
    'isdf_sigma_reduced_basis\([^;]+;\s*return', 'once')));
assert(exist(fullfile(repo_root, 'src', 'GW', 'ISDF', ...
    'isdf_sigma_reduced_basis.m'), 'file') == 0);
```

- [ ] **Step 2: Run the test and verify legacy-reference failure**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_sigma_shared_loop.m')"
```

Expected: FAIL because the legacy reduced workflow remains.

- [ ] **Step 3: Add reduced matrix-elements and contraction handles**

Replace the reduced error branch in `sigma_ops`:

```matlab
case 'reduced_basis'
    ops.matrix_elements = @(block) ...
        local_matrix_elements(ctx, block, true);
    ops.contract = @(block, matrix_elements) ...
        local_contract_reduced(ctx, block, matrix_elements);
```

`local_contract_reduced` moves `isdf_sigma_reduced_basis.m:168-240` into the strategy. For mapped screened-W:

```matlab
target_zeta = matrix_elements.space.zeta_g(1:block.n_cutoff, :);
kernel = isdf.screened_kernel( ...
    block.screened_w, target_zeta, block.coulg_cutoff);
for nn = 1:ctx.nbands
    coeff = matrix_elements.space.product_mu(:, nn);
    screened_value = ctx.fact * isdf.screened_contract(kernel, coeff);
    if block.occ_kq(nn) > 0
        asx_loc = asx_loc - block.occ_kq(nn) * screened_value;
    end
    ach_loc = ach_loc + screened_value;
end
```

For the full-inverse fallback, call the same static `sigma_cohsex` contraction as direct mode. Bare exchange always uses reconstructed `matrix_elements.gme`. For exact CH, construct the full screened matrix only when needed using `isdf.screened_kernel(screened_w, [], coulg_cutoff)`.

Change the reduced matrix-elements return value to:

```matlab
matrix_elements.gme = reshape( ...
    space.zeta_g * space.product_mu, numel(block.idx.q), ctx.nbands);
matrix_elements.space = space;
```

At the same time, update the Task 7 direct/matrix branches to return the same struct with `space=[]`:

```matlab
matrix_elements.gme = gme;
matrix_elements.space = [];
```

Update `local_contract_full` to read `matrix_elements.gme`. Update the shared loop's q=0 exact-CH assignment before contraction to:

```matlab
if ctx.sig.exact_static_ch && block.iq_fbz == 1
    aqsch{in, ispin} = matrix_elements.gme(:, in);
end
block.aqsch = aqsch;
```

This makes `matrix_elements` one stable struct type for all three modes; no contract implementation may rely on persistent or closure-mutated q state.

- [ ] **Step 4: Remove the reduced early dispatch and workflow file**

Delete `sigma.m:26-37` and let reduced mode enter the shared loop. Remove `src/GW/ISDF/isdf_sigma_reduced_basis.m` after its mapping/setup logic is present in `sigma_context` and its contraction logic is present in `sigma_ops`.

- [ ] **Step 5: Run reduced sigma regressions**

Run sequentially to avoid the known parallel MATLAB memory spike:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_sigma_shared_loop.m')"
matlab -wait -batch "run('test/ISDF/validation/test_AgBr_isdf_cauchy_cohsex_validation.m')"
matlab -wait -batch "run('test/ISDF/validation/test_AgBr_isdf_cauchy_cohsex_smw_validation.m')"
matlab -wait -batch "run('test/ISDF/validation/test_mos2_222_isdf_reduced_validation.m')"
```

Expected: all PASS; AgBr and MoS2 sigma/eqp0 relative errors remain below `1e-8`.

- [ ] **Step 6: Commit shared reduced sigma**

```powershell
git add -- src/GW/sigma/sigma.m src/GW/sigma/sigma_ops.m test/ISDF/unit/test_ISDF_sigma_shared_loop.m
git add -u -- src/GW/ISDF/isdf_sigma_reduced_basis.m
git diff --cached --check
git commit -m "Integrate reduced sigma into shared loop"
```

---

### Task 9: Move ISDF Implementations Behind the Package Boundary

**Files:**
- Modify: all public files under `src/GW/ISDF/+isdf/`
- Create: all locked files under `src/GW/ISDF/+isdf/private/`
- Modify: `src/GW/epsilon/epsilon_ops.m`
- Modify: `src/GW/sigma/sigma_ops.m`
- Modify: `test/ISDF/unit/test_ISDF_package_api.m`
- Modify: all `test/ISDF/unit/test_ISDF_*.m` callers of flat `isdf_*` numerical APIs
- Delete: superseded flat numerical files under `src/GW/ISDF/`

**Interfaces:**
- Consumes: package façade proven in Task 2.
- Produces: self-contained `+isdf` numerical implementation with private helpers and no repository call sites to superseded flat APIs.

- [ ] **Step 1: Add a failing package-boundary test**

Append to `test_ISDF_package_api.m`:

```matlab
flat_names = { ...
    'isdf_build_space', ...
    'isdf_reduced_polarizability', ...
    'isdf_static_screened_interaction', ...
    'isdf_screened_coulomb_kernel', ...
    'isdf_screened_coulomb_contract', ...
    'isdf_epsilon_batch', ...
    'isdf_sigma_batch', ...
    'isdf_wavefunction_real_component', ...
    'isdf_comega_cstar'};
for i = 1:numel(flat_names)
    assert(exist(flat_names{i}, 'file') == 0, ...
        'Superseded flat ISDF API remains: %s', flat_names{i});
end

assert(exist('isdf.build_space', 'file') == 2);
assert(exist('isdf.matrix_elements', 'file') == 2);
```

- [ ] **Step 2: Run the test and verify flat-API failure**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/unit/test_ISDF_package_api.m')"
```

Expected: FAIL because compatibility package files still call flat functions and the flat files exist.

- [ ] **Step 3: Move the product-space implementation into package/private files**

Implement package internals with these exact responsibilities:

- `build_space.m`: public scalar/component normalization, matrix-free vs explicit QRCP selection, construction of `ind_mu`, `product_mu`, `zeta_g`, selected components, rank and solve info.
- `private/component_products.m`: expose one callable private entry, not several inaccessible sibling subfunctions:

  ```matlab
  function [products, weight] = component_products( ...
      left, right, grid_indices, projection)
  ```

  Normalize scalar inputs to one-element cells. Empty `grid_indices` selects every real-space row; otherwise select exactly those rows. Empty `projection` returns explicit component products using the locked `conj(left).*right` convention; a nonempty projection returns the projected products used by matrix-free sampling. When requested, `weight` is the rowwise squared 2-norm summed over all component-pair columns. `build_space` is the only caller and chooses behavior through these data arguments—do not add an operation-string dispatcher.
- `private/product_gram.m`: move `isdf_product_gram` unchanged except for the shorter function name.
- `private/sample_points.m`: move QRCP, randomized QRCP, kmeans and weighted selection behind one `sample_points(left,right,options)` entry; retain existing RNG and weighting semantics.
- `private/stable_solve.m`: move `isdf_stable_right_solve` and warning restoration.
- `private/zeta_to_g.m`: merge `isdf_zeta_g_from_product_gram` and `isdf_zeta_real_to_g` into one private conversion path.

Do not change formulas or default values while moving code. Use the package test and existing component-space tests as equivalence gates after each file move.

- [ ] **Step 4: Move polarizability and screened-W implementations**

Replace façade bodies with the existing implementations:

- `polarizability.m` owns the current direct/Cauchy reduced coefficient code from `isdf_reduced_polarizability.m` and `isdf_comega_cstar.m` as local functions.
- `screened_w.m` owns the current SMW construction.
- `screened_kernel.m` owns the current kernel projection.
- `screened_contract.m` owns the current coefficient contraction.
- `real_component.m` owns the current FFT conversion.

Keep stable error identifiers and conjugation conventions. Update `matrix_elements.m` to call the now self-contained package `build_space`.

- [ ] **Step 5: Update repository call sites and unit tests**

Use:

```powershell
rg -n "isdf_(build_space|reduced_polarizability|static_screened_interaction|screened_coulomb_kernel|screened_coulomb_contract|epsilon_batch|sigma_batch|wavefunction_real_component|comega_cstar)" src test
```

Replace every source call with the locked package name. Rewrite numerical unit tests to test public package behavior rather than private helpers. For example:

```matlab
space = isdf.build_space(left_components, right_components, ...
    idx_q, fftgrid, options);
polar = isdf.polarizability(space, ev_occ, ev_unocc, solver_options);
```

No test may add `+isdf/private` to the MATLAB path.

- [ ] **Step 6: Delete superseded flat functions**

After `rg` returns no source/test callers, delete the superseded public workflow and numerical files, including the stale `isdf_kernelg_current_fft.asv`. Retain a flat helper only if `rg` shows a non-ISDF production caller; document that caller in the commit message.

- [ ] **Step 7: Run all ISDF unit tests**

Run every script returned by:

```powershell
Get-ChildItem test/ISDF/unit/test_ISDF_*.m | Select-Object -ExpandProperty FullName
```

Execute each with:

```powershell
matlab -wait -batch "run('<script-path>')"
```

Expected: every unit script exits with code 0.

- [ ] **Step 8: Verify the namespace boundary and commit**

Run:

```powershell
rg -n "isdf_(build_space|reduced_polarizability|static_screened_interaction|screened_coulomb_kernel|screened_coulomb_contract|epsilon_batch|sigma_batch|wavefunction_real_component|comega_cstar)" src test
git diff --check
```

Expected: `rg` prints no source/test matches and `git diff --check` exits 0.

Commit:

```powershell
git add -- src/GW/ISDF src/GW/epsilon/epsilon_ops.m src/GW/sigma/sigma_ops.m test/ISDF/unit
git diff --cached --check
git commit -m "Move ISDF kernels into MATLAB package"
```

---

### Task 10: Documentation, Full Regression, and Final Structural Verification

**Files:**
- Modify: `docs/isdf-interpolation-derivation.md`
- Modify: `docs/superpowers/specs/2026-07-24-isdf-shared-loop-refactor-design.md` only if final names differ from the locked plan because of a verified MATLAB restriction
- Test: all `test/ISDF/unit/*.m`
- Test: selected `test/ISDF/validation/*.m`

**Interfaces:**
- Consumes: completed shared-loop implementation.
- Produces: verified, documented refactor ready for code review.

- [ ] **Step 1: Update the derivation and usage documentation**

Document the final package calls:

```matlab
space = isdf.build_space(left_components, right_components, ...
    idx_q, fftgrid, isdf_options);
polar = isdf.polarizability(space, ev_occ, ev_unocc, solver_options);
screened = isdf.screened_w(space, coulg, polar);
```

Document that `epsilon.m` and `sigma.m` each own one shared loop and that direct/matrix/reduced differences enter through preselected ops.

- [ ] **Step 2: Run MATLAB Code Analyzer**

Run:

```powershell
matlab -wait -batch "files={'src/GW/epsilon/epsilon.m','src/GW/epsilon/epsilon_context.m','src/GW/epsilon/epsilon_prepare_block.m','src/GW/epsilon/epsilon_ops.m','src/GW/sigma/sigma.m','src/GW/sigma/sigma_context.m','src/GW/sigma/sigma_prepare_block.m','src/GW/sigma/sigma_ops.m'}; for k=1:numel(files), messages=checkcode(files{k},'-id'); assert(isempty(messages), files{k}); end"
```

Expected: exit code 0 with no analyzer messages.

- [ ] **Step 3: Run the full unit suite sequentially**

Run every `test/ISDF/unit/test_ISDF_*.m` script in a separate MATLAB process. Expected: all exit code 0.

- [ ] **Step 4: Run real-system validations sequentially**

Run:

```powershell
matlab -wait -batch "run('test/ISDF/validation/test_AgBr_isdf_epsilon_validation.m')"
matlab -wait -batch "run('test/ISDF/validation/test_AgBr_isdf_cauchy_epsilon_validation.m')"
matlab -wait -batch "run('test/ISDF/validation/test_AgBr_isdf_cauchy_cohsex_validation.m')"
matlab -wait -batch "run('test/ISDF/validation/test_AgBr_isdf_cauchy_cohsex_smw_validation.m')"
matlab -wait -batch "run('test/ISDF/validation/test_mos2_222_isdf_epsilon_validation.m')"
matlab -wait -batch "run('test/ISDF/validation/test_mos2_222_isdf_reduced_validation.m')"
```

Expected: every validation exits 0 without changing its existing tolerance. Run the two AgBr sigma validations sequentially, not in parallel.

- [ ] **Step 5: Verify structural acceptance criteria**

Run:

```powershell
rg -n "isdf_(epsilon|sigma)_reduced_basis|use_isdf|use_isdf_matrix_elements" src/GW/epsilon/epsilon.m src/GW/sigma/sigma.m
rg -n "isdf\.algorithm|isdf\.enable" src/GW/epsilon/epsilon.m src/GW/sigma/sigma.m
git diff --check
git status --short
```

Expected:

- the first two `rg` commands print no matches;
- `git diff --check` exits 0;
- status contains only the intended refactor/documentation files.

- [ ] **Step 6: Compare against the backup tag**

Run:

```powershell
git diff --stat isdf-before-shared-loop-refactor..HEAD
git show --no-patch --oneline isdf-before-shared-loop-refactor
```

Expected: the tag resolves to `e186aab`, and the diff shows shared-loop/package changes without unrelated directories.

- [ ] **Step 7: Commit documentation and verification updates**

```powershell
git add -- docs/isdf-interpolation-derivation.md docs/superpowers/specs/2026-07-24-isdf-shared-loop-refactor-design.md test/ISDF
git diff --cached --check
git commit -m "Document and verify shared ISDF workflows"
```

- [ ] **Step 8: Request final code review**

Use `superpowers:requesting-code-review` against the full diff from `isdf-before-shared-loop-refactor` to `HEAD`. The review must explicitly check shared-loop structure, low-rank memory preservation, scalar/spinor conjugation, symmetry mappings, exact CH, and public output compatibility.
