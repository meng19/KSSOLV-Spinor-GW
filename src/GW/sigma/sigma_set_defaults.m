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

if ~isfield(sig.isdf, 'enable')
    sig.isdf.enable = false;
end

if ~isfield(sig.isdf, 'rank_ratio')
    sig.isdf.rank_ratio = 1;
end

if ~isfield(sig.isdf, 'sample_method')
    sig.isdf.sample_method = 'qrcp';
end

if ~isfield(sig.isdf, 'seed')
    sig.isdf.seed = 0;
end

if ~isfield(sig.isdf, 'algorithm') || isempty(sig.isdf.algorithm)
    if sig.freq_dep == 0
        sig.isdf.algorithm = 'reduced_basis';
    else
        sig.isdf.algorithm = 'matrix_elements';
    end
end
if ~isfield(sig.isdf, 'reduced_solver') || isempty(sig.isdf.reduced_solver)
    sig.isdf.reduced_solver = 'cauchy';
end

if ~isfield(sig.isdf, 'cauchy_froErr') || isempty(sig.isdf.cauchy_froErr)
    sig.isdf.cauchy_froErr = 1e-8;
end

if ~isfield(sig.isdf, 'cauchy_MaxIter') || isempty(sig.isdf.cauchy_MaxIter)
    sig.isdf.cauchy_MaxIter = 12;
end
end
