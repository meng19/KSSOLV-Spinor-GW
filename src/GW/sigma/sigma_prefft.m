function idx = sigma_prefft(wfnkq, wfnk, q_mtx, iq, ik, sys, idx, use_gpu)
% For iq, precompute indices of wav in fft_grid in sigma with GPU support
% 支持两种模式：
%   1. 预计算模式（precompute_wav=true）：idx 是 cell 数组
%   2. 按需计算模式（precompute_wav=false）：idx 是单个 struct（节省内存）

%% set up fft_grid
Nfft = [sys.n1, sys.n2, sys.n3];
grid_size = [Nfft(1), Nfft(2), Nfft(3)];

%% 判断当前模式：是否为按需计算（on-the-fly）
is_on_the_fly = isempty(idx) || ~isstruct(idx) || ~isfield(idx, 'kq') || ~iscell(idx.kq);

%% gvec to fft grid
g_kq = wfnkq.mill;
g_k  = wfnk.mill;
g_q  = q_mtx;

bidx_kq = g2fft_grid(g_kq, Nfft(1), Nfft(2), Nfft(3));
bidx_k  = g2fft_grid(g_k,  Nfft(1), Nfft(2), Nfft(3));
bidx_q  = g2fft_grid(g_q,  Nfft(1), Nfft(2), Nfft(3));

%% Compute linear indices
lin_idx_kq = bidx_kq(:,1) + (bidx_kq(:,2)-1)*grid_size(1) ...
           + (bidx_kq(:,3)-1)*grid_size(1)*grid_size(2);
       
lin_idx_k  = bidx_k(:,1)  + (bidx_k(:,2)-1)*grid_size(1)  ...
           + (bidx_k(:,3)-1)*grid_size(1)*grid_size(2);
       
lin_idx_q  = bidx_q(:,1)  + (bidx_q(:,2)-1)*grid_size(1)  ...
           + (bidx_q(:,3)-1)*grid_size(1)*grid_size(2);

%% Store indices according to mode
if is_on_the_fly
    % ==================== 按需计算模式（precompute_wav=false） ====================
    % 返回单个 struct，而不是 cell
    idx.k  = lin_idx_k;
    idx.q  = lin_idx_q;
    idx.kq = lin_idx_kq;
    
    if use_gpu
        try
            idx.k  = gpuArray(idx.k);
            idx.q  = gpuArray(idx.q);
            idx.kq = gpuArray(idx.kq);
        catch
            warning('Failed to move indices to GPU, using CPU.');
        end
    end
    
else
    % ==================== 预计算模式（precompute_wav=true） ====================
    if use_gpu
        try
            idx.kq{iq, ik} = gpuArray(lin_idx_kq);
            idx.k{ik}      = gpuArray(lin_idx_k);
            idx.q{iq, ik}  = gpuArray(lin_idx_q);
        catch ME
            if strcmp(ME.identifier, 'parallel:gpu:array:OOM')
                warning('GPU memory insufficient for sigma indices. Falling back to CPU.');
                idx.kq{iq, ik} = lin_idx_kq;
                idx.k{ik}      = lin_idx_k;
                idx.q{iq, ik}  = lin_idx_q;
            else
                rethrow(ME);
            end
        end
    else
        % CPU version
        idx.kq{iq, ik} = lin_idx_kq;
        idx.k{ik}      = lin_idx_k;
        idx.q{iq, ik}  = lin_idx_q;
    end
end

%% Store grid information (只存一次)
if is_on_the_fly
    if ~isfield(idx, 'grid_size')
        idx.grid_size = grid_size;
    end
    if ~isfield(idx, 'use_gpu')
        idx.use_gpu = use_gpu && isa(idx.kq, 'gpuArray');
    end
else
    if ~isfield(idx, 'grid_size')
        idx.grid_size = grid_size;
    end
    if ~isfield(idx, 'use_gpu')
        idx.use_gpu = use_gpu && isa(idx.kq{1}, 'gpuArray');  % 检查第一个元素
    end
end

end