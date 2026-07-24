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

if ~isfield(eps.isdf, 'enable')
    eps.isdf.enable = false;
end

if ~isfield(eps.isdf, 'algorithm') || isempty(eps.isdf.algorithm)
    if eps.freq_dep == 0
        eps.isdf.algorithm = 'reduced_basis';
    else
        eps.isdf.algorithm = 'matrix_elements';
    end
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

if ~isfield(eps.isdf, 'reduced_solver') || isempty(eps.isdf.reduced_solver)
    eps.isdf.reduced_solver = 'cauchy';
end

if ~isfield(eps.isdf, 'output') || isempty(eps.isdf.output)
    eps.isdf.output = 'screened_w';
end

if ~isfield(eps.isdf, 'cauchy_froErr') || isempty(eps.isdf.cauchy_froErr)
    eps.isdf.cauchy_froErr = 1e-8;
end

if ~isfield(eps.isdf, 'cauchy_MaxIter') || isempty(eps.isdf.cauchy_MaxIter)
    eps.isdf.cauchy_MaxIter = 12;
end

end
