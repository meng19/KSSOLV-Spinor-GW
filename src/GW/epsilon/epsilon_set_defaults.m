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

if ~isfield(eps, 'isdf') || isempty(eps.isdf)
    eps.isdf.enable = false;
end

if ~isfield(eps.isdf, 'enable')
    eps.isdf.enable = false;
end

if ~isfield(eps.isdf, 'rank_ratio')
    eps.isdf.rank_ratio = 1;
end

if ~isfield(eps.isdf, 'sample_method')
    eps.isdf.sample_method = 'qrcp';
end

if ~isfield(eps.isdf, 'seed')
    eps.isdf.seed = 0;
end
end
