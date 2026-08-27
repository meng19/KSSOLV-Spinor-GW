script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

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
options.warn_rank_selection = false;

space = isdf.build_space(conj(phi), psi, idx_q, fftgrid, options);

assert(isfield(space, 'ind_mu'));
assert(isfield(space, 'zeta_g'));
assert(isfield(space, 'phi_mu'));
assert(isfield(space, 'psi_mu'));
assert(numel(space.ind_mu) == options.rank);
assert(size(space.zeta_g, 2) == options.rank);
assert(size(space.phi_mu, 1) == options.rank);
assert(size(space.psi_mu, 1) == options.rank);

auto_options = rmfield(options, {'rank', 'warn_rank_selection'});
auto_options.rank_ratio = 1;
auto_space = isdf.build_space(conj(phi), psi, idx_q, fftgrid, ...
    auto_options);
expected_rank = ceil(sqrt(nphi * npsi));
assert(auto_space.rank == expected_rank);
assert(auto_space.options.recommended_rank == expected_rank);
assert(strcmp(auto_space.options.rank_source, 'auto'));

warning_state = warning('query', 'ISDF:RankAboveRecommended');
warning_cleanup = onCleanup(@() warning( ...
    warning_state.state, 'ISDF:RankAboveRecommended'));
warning('on', 'ISDF:RankAboveRecommended');
rank_warning_options = options;
rank_warning_options.warn_rank_selection = true;
rank_warning_options.rank_ratio = 0.987;
lastwarn('');
isdf.build_space(conj(phi), psi, idx_q, fftgrid, rank_warning_options);
[~, warning_id] = lastwarn;
assert(strcmp(warning_id, 'ISDF:RankAboveRecommended'), ...
    'Explicit rank above the ISDF default did not emit a warning.');

low_warning_state = warning('query', 'ISDF:RankBelowRecommended');
low_warning_cleanup = onCleanup(@() warning( ...
    low_warning_state.state, 'ISDF:RankBelowRecommended'));
warning('on', 'ISDF:RankBelowRecommended');
low_rank_options = options;
low_rank_options.rank = 1;
low_rank_options.warn_rank_selection = true;
lastwarn('');
isdf.build_space(conj(phi), psi, idx_q, fftgrid, low_rank_options);
[~, warning_id] = lastwarn;
assert(strcmp(warning_id, 'ISDF:RankBelowRecommended'), ...
    'Explicit rank below the ISDF default did not emit a warning.');

fprintf('ISDF build space test passed.\n');
