function eps = epsilon_set_defaults(eps)
% ==================== 设置 sig 的默认值 ====================
if ~isfield(eps, 'use_gpu')
    eps.use_gpu = false;                % 默认不使用GPU
end

if ~isfield(eps, 'save_mem')
    eps.save_mem = false;
end

if ~isfield(eps, 'precompute_wav')
    eps.precompute_wav = false;
end

% A non-empty checkpoint file saves the completed screening data together
% with the inputs needed to resume at sigma in a later MATLAB session.
if ~isfield(eps, 'checkpoint_file') || isempty(eps.checkpoint_file)
    eps.checkpoint_file = '';
end
if isfield(eps, 'checkpoint') && ~isempty(eps.checkpoint)
    eps.checkpoint_file = eps.checkpoint;
elseif isfield(eps, 'breakpoint') && ~isempty(eps.breakpoint)
    eps.checkpoint_file = eps.breakpoint;
end
if ~(ischar(eps.checkpoint_file) || ...
        (isstring(eps.checkpoint_file) && isscalar(eps.checkpoint_file)))
    error('GW:EpsilonCheckpointFile', ...
        'eps.checkpoint_file must be a character vector or string scalar.');
end

if ~isfield(eps, 'freq_dep')
    eps.freq_dep = 0;                  % 默认静态COHSEX
end

if ~isfield(eps, 'freq_dep_method')
    eps.freq_dep_method = 2;           % 默认方法
end

if ~isfield(eps, 'isdf') || isempty(eps.isdf) || ~isstruct(eps.isdf)
    eps.isdf = struct();
end
eps.isdf = isdf.workflow_defaults(eps.isdf, eps.freq_dep);

if ~isfield(eps.isdf, 'output') || isempty(eps.isdf.output)
    eps.isdf.output = 'screened_w';
end
if ~isfield(eps.isdf, 'screened_w_storage') || ...
        isempty(eps.isdf.screened_w_storage)
    % Sigma uses zeta_g, epsilon_vcoul, and k_mu only.  Do not retain the
    % O(rank^2 * nfreq) polar diagnostic pages in production by default.
    eps.isdf.screened_w_storage = 'minimal';
end
eps.isdf.screened_w_storage = lower(char(eps.isdf.screened_w_storage));
if ~any(strcmp(eps.isdf.screened_w_storage, {'minimal', 'full'}))
    error('ISDF:ScreenedWStorage', ...
        ['eps.isdf.screened_w_storage must be ''minimal'' or ''full'' ' ...
         '(got ''%s'').'], eps.isdf.screened_w_storage);
end
if ~isfield(eps.isdf, 'cache_real_wfn') || isempty(eps.isdf.cache_real_wfn)
    eps.isdf.cache_real_wfn = false;
end

end
