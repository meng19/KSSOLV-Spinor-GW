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
if ~isfield(eps.isdf, 'cache_real_wfn') || isempty(eps.isdf.cache_real_wfn)
    eps.isdf.cache_real_wfn = false;
end

end
