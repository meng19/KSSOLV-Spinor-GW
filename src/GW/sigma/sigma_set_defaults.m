function sig = sigma_set_defaults(sig)
% ==================== 设置 sig 的默认值 ====================
if ~isfield(sig, 'use_gpu')
    sig.use_gpu = false;                % 默认不使用GPU
end

if ~isfield(sig, 'precompute_wav')
    sig.precompute_wav = false;
end

if ~isfield(sig, 'freq_dep')
    sig.freq_dep = 0;                  % 默认静态COHSEX
end

if ~isfield(sig, 'freq_dep_method')
    sig.freq_dep_method = 2;           % 默认方法
end

if ~isfield(sig, 'no_symmetries_q_grid')
    sig.no_symmetries_q_grid = false;
end

if ~isfield(sig, 'exact_static_ch')
    sig.exact_static_ch = false;
end

if ~isfield(sig, 'isdf') || isempty(sig.isdf) || ~isstruct(sig.isdf)
    sig.isdf = struct();
end
sig.isdf = isdf.workflow_defaults(sig.isdf, sig.freq_dep);

% The legacy path uses the NN product space for every sigma contribution.
% A VN space can be selected explicitly for the bare-exchange term only.
if ~isfield(sig.isdf, 'exchange_space') || isempty(sig.isdf.exchange_space)
    sig.isdf.exchange_space = 'nn';
end
if ~isfield(sig.isdf, 'global_nn_space') || isempty(sig.isdf.global_nn_space)
    sig.isdf.global_nn_space = false;
end
if ~isfield(sig.isdf, 'global_vn_space') || isempty(sig.isdf.global_vn_space)
    sig.isdf.global_vn_space = false;
end
if ~isfield(sig.isdf, 'reuse_nn_for_vn') || isempty(sig.isdf.reuse_nn_for_vn)
    sig.isdf.reuse_nn_for_vn = false;
end
if ~isfield(sig.isdf, 'reuse_eps_real_wfn') || isempty(sig.isdf.reuse_eps_real_wfn)
    sig.isdf.reuse_eps_real_wfn = false;
end

if ~isfield(sig.isdf, 'screened_kernel_cache_gb') || ...
        isempty(sig.isdf.screened_kernel_cache_gb)
    % Cache budget for reusing reduced screened kernels shared by several
    % diagonal bands.  Inf keeps every reusable kernel (no memory cap); a
    % non-positive value disables the reuse.  The value is in GB.
    if isfield(sig.isdf, 'screened_kernel_cache_bytes') && ...
            ~isempty(sig.isdf.screened_kernel_cache_bytes)
        sig.isdf.screened_kernel_cache_gb = ...
            sig.isdf.screened_kernel_cache_bytes / (1024^3);
    else
        sig.isdf.screened_kernel_cache_gb = Inf;
    end
end
if isfield(sig.isdf, 'screened_kernel_cache_bytes')
    % Backward compatibility: keep legacy alias in bytes, but prefer the GB
    % field for user input.
    sig.isdf.screened_kernel_cache_bytes = ...
        sig.isdf.screened_kernel_cache_gb * (1024^3);
end
if ~(isnumeric(sig.isdf.screened_kernel_cache_gb) && ...
        isscalar(sig.isdf.screened_kernel_cache_gb) && ...
        sig.isdf.screened_kernel_cache_gb >= 0)
    error('ISDF:ScreenedKernelCacheLimit', ...
        ['sig.isdf.screened_kernel_cache_gb must be a ' ...
         'non-negative scalar or Inf.']);
end
end
