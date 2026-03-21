function sig = sigma_set_defaults(sig)
% ==================== 设置 sig 的默认值 ====================
if ~isfield(sig, 'use_gpu')
    sig.use_gpu = false;                % 默认不使用GPU
end

if ~isfield(sig, 'freq_dep')
    sig.freq_dep = 0;                  % 默认静态COHSEX
end

if ~isfield(sig, 'freq_dep_method')
    sig.freq_dep_method = 2;           % 默认方法
end
end