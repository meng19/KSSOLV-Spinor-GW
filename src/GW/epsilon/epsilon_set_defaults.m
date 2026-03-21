function eps = sigma_set_defaults(eps)
% ==================== 设置 sig 的默认值 ====================
if ~isfield(eps, 'use_gpu')
    eps.use_gpu = false;                % 默认不使用GPU
end

if ~isfield(eps, 'freq_dep')
    eps.freq_dep = 0;                  % 默认静态COHSEX
end

if ~isfield(eps, 'freq_dep_method')
    eps.freq_dep_method = 2;           % 默认方法
end
end