# ISDF Cauchy COHSEX Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an optional static COHSEX path that uses ISDF spaces plus Cauchy-integral decoupling to reduce the dominant polarizability/self-energy contractions.

**Architecture:** Keep the existing `epsilon` and `sigma` defaults unchanged. Add focused helpers under `src/GW/ISDF/`, then dispatch from `sigma.m` only when `sig.isdf.enable = true` and `sig.isdf.algorithm = 'cauchy_cohsex'`. Validate the first implementation against existing direct static COHSEX on AgBr/Gamma with full ISDF rank.

**Tech Stack:** MATLAB, existing KSSOLV GW data structures, existing `src/GW/ISDF` helpers, existing MATLAB validation scripts in `test/ISDF`.

## Global Constraints

- Do not change default `epsilon` or `sigma` behavior.
- First implementation supports `eps.freq_dep = 0` and `sig.freq_dep = 0` only.
- First implementation supports CPU only.
- First implementation supports Gamma/single-k validation first.
- Unsupported modes must fail explicitly.
- Tests compare against the existing direct path because the user considers current calculations correct.

---

### Task 1: Add Cauchy COHSEX Defaults And Dispatch Guard

**Files:**
- Modify: `src/GW/sigma/sigma_set_defaults.m`
- Modify: `src/GW/sigma/sigma.m`
- Test: `test/ISDF/test_ISDF_cauchy_cohsex_guards.m`

**Interfaces:**
- Consumes: existing `sig.isdf` struct.
- Produces: `sig.isdf.algorithm`, `sig.isdf.cauchy_method`, `sig.isdf.cauchy_froErr`, `sig.isdf.cauchy_MaxIter`; dispatch to `isdf_sigma_cohsex_cauchy`.

- [ ] **Step 1: Add a failing guard test**

Create `test/ISDF/test_ISDF_cauchy_cohsex_guards.m`:

```matlab
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

sig = struct();
sig = sigma_set_defaults(sig);
assert(isfield(sig, 'isdf'));
assert(isfield(sig.isdf, 'algorithm'));
assert(strcmp(sig.isdf.algorithm, 'matrix_elements'));
assert(isfield(sig.isdf, 'cauchy_method'));
assert(strcmp(sig.isdf.cauchy_method, 'cauchy'));
assert(isfield(sig.isdf, 'cauchy_froErr'));
assert(isfield(sig.isdf, 'cauchy_MaxIter'));

fprintf('ISDF Cauchy COHSEX default guard test passed.\n');
```

- [ ] **Step 2: Run the guard test and verify it fails**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_cauchy_cohsex_guards.m')"""
```

Expected: FAIL because the new default fields do not exist yet.

- [ ] **Step 3: Add sigma defaults**

In `src/GW/sigma/sigma_set_defaults.m`, after the existing `sig.isdf` defaults:

```matlab
if ~isfield(sig.isdf, 'algorithm') || isempty(sig.isdf.algorithm)
    sig.isdf.algorithm = 'matrix_elements';
end

if ~isfield(sig.isdf, 'cauchy_method') || isempty(sig.isdf.cauchy_method)
    sig.isdf.cauchy_method = 'cauchy';
end

if ~isfield(sig.isdf, 'cauchy_froErr') || isempty(sig.isdf.cauchy_froErr)
    sig.isdf.cauchy_froErr = 1e-8;
end

if ~isfield(sig.isdf, 'cauchy_MaxIter') || isempty(sig.isdf.cauchy_MaxIter)
    sig.isdf.cauchy_MaxIter = 12;
end
```

- [ ] **Step 4: Add dispatch guard in `sigma.m`**

Near the beginning of `src/GW/sigma/sigma.m`, after defaults and GPU/use_isdf setup:

```matlab
if use_isdf && isfield(sig.isdf, 'algorithm') && strcmpi(sig.isdf.algorithm, 'cauchy_cohsex')
    if use_gpu
        error('ISDF:CauchyCOHSEXGPU', ...
            'ISDF Cauchy COHSEX path currently supports CPU execution only.');
    end
    if sig.freq_dep ~= 0
        error('ISDF:CauchyCOHSEXFrequency', ...
            'ISDF Cauchy COHSEX path requires sig.freq_dep = 0.');
    end
    sig = isdf_sigma_cohsex_cauchy(eps, sig, sys, options, syms);
    return;
end
```

- [ ] **Step 5: Run guard test and existing ISDF matrix tests**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_cauchy_cohsex_guards.m')"""
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_comega_cstar.m')"""
```

Expected: PASS for both.

---

### Task 2: Build A Reusable ISDF Space Object

**Files:**
- Create: `src/GW/ISDF/isdf_build_space.m`
- Test: `test/ISDF/test_ISDF_build_space.m`

**Interfaces:**
- Consumes: real-space matrices `phi`, `psi`, FFT index `idx_q`, FFT grid, and ISDF options.
- Produces: struct with fields `phi`, `psi`, `ind_mu`, `zeta_g`, `phi_mu`, `psi_mu`, `rank`.

- [ ] **Step 1: Add a failing unit test**

Create `test/ISDF/test_ISDF_build_space.m`:

```matlab
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(repo_root);

rng(11);
ngrid = 16;
nphi = 3;
npsi = 4;
phi = randn(ngrid, nphi) + 1i * randn(ngrid, nphi);
psi = randn(ngrid, npsi) + 1i * randn(ngrid, npsi);
idx_q = (1:ngrid).';
fftgrid = [4, 4, 1];

options = struct();
options.rank = nphi * npsi;
options.sample_method = 'qrcp';
options.seed = 0;

space = isdf_build_space(phi, psi, idx_q, fftgrid, options);

assert(isfield(space, 'ind_mu'));
assert(isfield(space, 'zeta_g'));
assert(isfield(space, 'phi_mu'));
assert(isfield(space, 'psi_mu'));
assert(numel(space.ind_mu) == options.rank);
assert(size(space.zeta_g, 2) == options.rank);
assert(size(space.phi_mu, 1) == options.rank);
assert(size(space.psi_mu, 1) == options.rank);

fprintf('ISDF build space test passed.\n');
```

- [ ] **Step 2: Run the test and verify it fails**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_build_space.m')"""
```

Expected: FAIL because `isdf_build_space` does not exist.

- [ ] **Step 3: Implement `isdf_build_space.m`**

Create `src/GW/ISDF/isdf_build_space.m`:

```matlab
function space = isdf_build_space(phi, psi, idx_q, fftgrid, options)
%ISDF_BUILD_SPACE Build a compact ISDF product-space representation.

if nargin < 5 || isempty(options)
    options = struct();
end

options = isdf_set_defaults(options, size(phi, 2), size(psi, 2), size(phi, 1));
ind_mu = isdf_indices(phi, psi, options);
zeta_g = isdf_kernelg_current_fft(phi, psi, ind_mu, idx_q, fftgrid);

space = struct();
space.phi = phi;
space.psi = psi;
space.ind_mu = ind_mu;
space.zeta_g = zeta_g;
space.phi_mu = phi(ind_mu, :);
space.psi_mu = psi(ind_mu, :);
space.rank = numel(ind_mu);
space.options = options;
end
```

- [ ] **Step 4: Run the test**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_build_space.m')"""
```

Expected: PASS.

---

### Task 3: Wrap Cauchy Polarizability In Reduced Space

**Files:**
- Create: `src/GW/ISDF/isdf_cauchy_polarizability.m`
- Modify: `src/GW/ISDF/isdf_comega_cstar.m`
- Test: `test/ISDF/test_ISDF_cauchy_polarizability.m`

**Interfaces:**
- Consumes: `vc_space`, occupied energies, conduction energies, Cauchy options.
- Produces: `polar.coeff`, `polar.info`, where `coeff = C * Omega^{-1} * C'`.

- [ ] **Step 1: Add a failing unit test**

Create `test/ISDF/test_ISDF_cauchy_polarizability.m`:

```matlab
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(repo_root);

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
direct = isdf_cauchy_polarizability(vc_space, ev_occ, ev_unocc, direct_options);

cauchy_options = struct('method', 'cauchy', 'froErr', 1e-10, 'MaxIter', 8);
cauchy = isdf_cauchy_polarizability(vc_space, ev_occ, ev_unocc, cauchy_options);

relerr = norm(direct.coeff - cauchy.coeff, 'fro') / max(1, norm(direct.coeff, 'fro'));
assert(relerr < 1e-8, 'Cauchy polarizability differs from direct: %.3e', relerr);

fprintf('ISDF Cauchy polarizability test passed. relerr = %.3e\n', relerr);
```

- [ ] **Step 2: Run the test and verify it fails**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_cauchy_polarizability.m')"""
```

Expected: FAIL because `isdf_cauchy_polarizability` does not exist.

- [ ] **Step 3: Implement `isdf_cauchy_polarizability.m`**

Create `src/GW/ISDF/isdf_cauchy_polarizability.m`:

```matlab
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
```

- [ ] **Step 4: Remove dead elliptic helper functions from `isdf_comega_cstar.m`**

Delete unused local functions below `isdf_comega_cstar_cauchy`:

- `isdf_cauchy_integrand`
- `isdf_ellipk_modulus`
- `isdf_ellipjc`
- `isdf_ellipkkp`

The current implementation uses a circular contour and does not call these helpers.

- [ ] **Step 5: Run tests**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_cauchy_polarizability.m')"""
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_comega_cstar.m')"""
```

Expected: PASS.

---

### Task 4: Add Static Reduced Screened Interaction Helper

**Files:**
- Create: `src/GW/ISDF/isdf_static_screened_interaction.m`
- Test: `test/ISDF/test_ISDF_static_screened_interaction.m`

**Interfaces:**
- Consumes: `vc_space.zeta_g`, Coulomb vector `vcoul`, and `polar.coeff`.
- Produces: struct with reduced matrices `eps_mu`, `eps_mu_inv`, `w_mu`.

- [ ] **Step 1: Add a failing algebraic test**

Create `test/ISDF/test_ISDF_static_screened_interaction.m`:

```matlab
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(repo_root);

rng(13);
ng = 8;
nmu = 4;
vc_space = struct();
vc_space.zeta_g = randn(ng, nmu) + 1i * randn(ng, nmu);
vcoul = abs(randn(ng, 1)) + 0.5;
a = randn(nmu);
polar = struct();
polar.coeff = -(a * a' + eye(nmu));

screened = isdf_static_screened_interaction(vc_space, vcoul, polar);

vmat = vc_space.zeta_g' * diag(vcoul) * vc_space.zeta_g;
expected_eps_mu = inv(polar.coeff) - vmat;
expected_eps_mu_inv = inv(expected_eps_mu);
expected_w_mu = expected_eps_mu_inv;

assert(norm(screened.eps_mu - expected_eps_mu, 'fro') < 1e-10);
assert(norm(screened.eps_mu_inv - expected_eps_mu_inv, 'fro') < 1e-10);
assert(norm(screened.w_mu - expected_w_mu, 'fro') < 1e-10);

fprintf('ISDF static screened interaction test passed.\n');
```

- [ ] **Step 2: Run the test and verify it fails**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_static_screened_interaction.m')"""
```

Expected: FAIL because `isdf_static_screened_interaction` does not exist.

- [ ] **Step 3: Implement `isdf_static_screened_interaction.m`**

Create `src/GW/ISDF/isdf_static_screened_interaction.m`:

```matlab
function screened = isdf_static_screened_interaction(vc_space, vcoul, polar)
%ISDF_STATIC_SCREENED_INTERACTION Construct reduced static screened operator.

vcoul = vcoul(:);
if size(vc_space.zeta_g, 1) ~= numel(vcoul)
    error('ISDF:ScreenedInteractionSize', ...
        'Coulomb vector length must match zeta_g row count.');
end

vmat = vc_space.zeta_g' * (vcoul .* vc_space.zeta_g);
eps_mu = inv(polar.coeff) - vmat;
eps_mu_inv = inv(eps_mu);

screened = struct();
screened.vmat = vmat;
screened.eps_mu = eps_mu;
screened.eps_mu_inv = eps_mu_inv;
screened.w_mu = eps_mu_inv;
end
```

- [ ] **Step 4: Run the test**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_static_screened_interaction.m')"""
```

Expected: PASS.

---

### Task 5: Implement Static Cauchy COHSEX Sigma Path

**Files:**
- Create: `src/GW/ISDF/isdf_sigma_cohsex_cauchy.m`
- Modify: `src/GW/sigma/sigma.m`
- Test: `test/ISDF/test_AgBr_isdf_cauchy_cohsex_validation.m`

**Interfaces:**
- Consumes: existing `eps`, `sig`, `sys`, `options`, `syms`.
- Produces: a sigma-compatible output struct with `sig.sig` and `sig.eqp0`.

- [ ] **Step 1: Add AgBr validation test**

Create `test/ISDF/test_AgBr_isdf_cauchy_cohsex_validation.m` by adapting the existing AgBr static sigma validation:

```matlab
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

run('test_AgBr_spinor_gw.m');

eps_input = epsilon_set_defaults(eps);
eps_input.freq_dep = 0;
eps_result = epsilon(sys, options, syms, eps_input);

sig_base = sigma_set_defaults(sig);
sig_base.freq_dep = 0;

sig_direct = sigma(eps_result, sig_base, sys, options, syms);

sig_cauchy_input = sig_base;
sig_cauchy_input.isdf.enable = true;
sig_cauchy_input.isdf.algorithm = 'cauchy_cohsex';
sig_cauchy_input.isdf.sample_method = 'qrcp';
sig_cauchy_input.isdf.rank = sig_base.nbnd;
sig_cauchy_input.isdf.cauchy_method = 'cauchy';
sig_cauchy_input.isdf.cauchy_froErr = 1e-8;
sig_cauchy_input.isdf.cauchy_MaxIter = 12;

sig_cauchy = sigma(eps_result, sig_cauchy_input, sys, options, syms);

sig_relative_error = norm(sig_direct.sig(:) - sig_cauchy.sig(:)) / max(1, norm(sig_direct.sig(:)));
eqp_relative_error = norm(sig_direct.eqp0(:) - sig_cauchy.eqp0(:)) / max(1, norm(sig_direct.eqp0(:)));

fprintf('AgBr Cauchy COHSEX sig relative error = %.3e\n', sig_relative_error);
fprintf('AgBr Cauchy COHSEX eqp0 relative error = %.3e\n', eqp_relative_error);

assert(sig_relative_error < 1e-6, ...
    'AgBr Cauchy COHSEX sig validation failed: %.3e', sig_relative_error);
assert(eqp_relative_error < 1e-6, ...
    'AgBr Cauchy COHSEX eqp0 validation failed: %.3e', eqp_relative_error);
```

- [ ] **Step 2: Run the validation test and verify it fails**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_AgBr_isdf_cauchy_cohsex_validation.m')"""
```

Expected: FAIL because `isdf_sigma_cohsex_cauchy` is not implemented.

- [ ] **Step 3: Implement a guarded skeleton**

Create `src/GW/ISDF/isdf_sigma_cohsex_cauchy.m`:

```matlab
function sig = isdf_sigma_cohsex_cauchy(eps, sig, sys, options, syms)
%ISDF_SIGMA_COHSEX_CAUCHY Static COHSEX self-energy with ISDF/Cauchy.

if eps.freq_dep ~= 0 || sig.freq_dep ~= 0
    error('ISDF:CauchyCOHSEXFrequency', ...
        'ISDF Cauchy COHSEX requires eps.freq_dep = 0 and sig.freq_dep = 0.');
end

if isfield(sig, 'use_gpu') && sig.use_gpu
    error('ISDF:CauchyCOHSEXGPU', ...
        'ISDF Cauchy COHSEX path currently supports CPU execution only.');
end

if sys.nkpts ~= 1
    error('ISDF:CauchyCOHSEXKPoints', ...
        'Initial ISDF Cauchy COHSEX implementation supports single-k systems only.');
end

sig = isdf_sigma_cohsex_cauchy_singlek(eps, sig, sys, options, syms);
end
```

In the same file, add a private `isdf_sigma_cohsex_cauchy_singlek` function that initially calls the existing direct `sigma` path with `sig.isdf.algorithm = 'matrix_elements'` removed only as a temporary scaffold. The next step replaces the scaffold with the reduced formula. This scaffold must be removed before Task 5 is considered complete.

- [ ] **Step 4: Replace scaffold with reduced static formula**

Implement `isdf_sigma_cohsex_cauchy_singlek` using these operations:

```matlab
% 1. Load Gamma k wavefunctions using existing GW loaders.
% 2. Build real-space occupied, conduction, target, and summation wavefunction matrices.
% 3. Build vc_space = isdf_build_space(conj(occ_real), cond_real, idx.q, fftgrid, vc_options).
% 4. Build vn_space = isdf_build_space(conj(occ_real), target_real, idx.q, fftgrid, vn_options).
% 5. Build nn_space = isdf_build_space(conj(sum_real), target_real, idx.q, fftgrid, nn_options).
% 6. polar = isdf_cauchy_polarizability(vc_space, ev_occ, ev_cond, cauchy_options).
% 7. screened = isdf_static_screened_interaction(vc_space, vcoul, polar).
% 8. Contract reduced W with vn_space and nn_space to form SEX and COH.
% 9. Fill sig.sig and sig.eqp0 in the same shape as existing sigma.m output.
```

The contraction should follow the thesis Chapter 3 reduced-space form and the existing reference implementation from the old GW submodule:

```matlab
epsvc = vc_space.zeta_g / screened.eps_mu;
epsvcDcoul_vn = epsvc' * (vcoul .* vn_space.zeta_g);
epsvcDcoul_nn = epsvc' * (vcoul .* nn_space.zeta_g);
vcDcoul_vn = vc_space.zeta_g' * (vcoul .* vn_space.zeta_g);
vcDcoul_nn = vc_space.zeta_g' * (vcoul .* nn_space.zeta_g);

w1_vn = vcDcoul_vn' * epsvcDcoul_vn;
w1_nn = vcDcoul_nn' * epsvcDcoul_nn;
```

Use existing direct `sigma.m` normalization and band-energy update formulas as the reference for signs, volume factors, and `eqp0`.

- [ ] **Step 5: Run AgBr validation**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_AgBr_isdf_cauchy_cohsex_validation.m')"""
```

Expected: PASS with `sig_relative_error < 1e-6` and `eqp_relative_error < 1e-6`.

---

### Task 6: Regression Test Existing ISDF Paths

**Files:**
- No new files.
- Test: existing `test/ISDF/*.m` validation scripts.

**Interfaces:**
- Consumes: all previous implementation tasks.
- Produces: evidence that existing matrix-element ISDF and the new Cauchy path coexist.

- [ ] **Step 1: Run unit-level ISDF tests**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_epsilon_matrix_elements.m')"""
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_sigma_matrix_elements.m')"""
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_comega_cstar.m')"""
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_build_space.m')"""
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_cauchy_polarizability.m')"""
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_ISDF_static_screened_interaction.m')"""
```

Expected: PASS for all.

- [ ] **Step 2: Run AgBr validation tests**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_AgBr_isdf_epsilon_validation.m')"""
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_AgBr_isdf_sigma_validation.m')"""
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_AgBr_isdf_cauchy_cohsex_validation.m')"""
```

Expected: PASS for all.

- [ ] **Step 3: Run MoS2 existing ISDF validation**

Run:

```powershell
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_mos2_222_isdf_epsilon_validation.m')"""
C:\WINDOWS\system32\cmd.exe /c "matlab -wait -batch ""run('test/ISDF/test_mos2_222_isdf_sigma_validation.m')"""
```

Expected: PASS for existing matrix-element ISDF path. The new `cauchy_cohsex` path is not required to support MoS2 multi-k in this plan.

---

## Self-Review

- Spec coverage: the plan adds explicit defaults, a dispatch path, reduced ISDF space construction, Cauchy polarizability, static reduced screened interaction, and an AgBr validation against the existing direct path.
- Scope check: full-frequency GW/CD, GPU execution, and multi-k Cauchy COHSEX are excluded from this implementation.
- Type consistency: all new public helper functions use MATLAB structs and existing ISDF helper conventions.
- Testability: each new helper has a unit test, and the final path has an AgBr end-to-end validation.
