function [fft, idx] = epsilon_prefft(wfnkq, wfnk, iq, ik, pol, fft, idx, use_gpu)
% For iq, precompute indices of wav in fft_grid in epsilon with GPU support
% 支持两种模式：
%   1. 预计算模式（fft 和 idx 是 cell 数组）—— 用于 save_mem=false
%   2. 按需计算模式（fft 和 idx 是单个 struct）—— 用于 save_mem=true

%% set up fft_grid
Nfac = 3;
Nfft = zeros(1, 3);
for i = 1:3
    Nfft(i) = pol.fftgrid{iq}(1, i);
    while ~check_FFT_size(Nfft(i), Nfac)
        Nfft(i) = Nfft(i) + 1;
    end
end

grid_size = [Nfft(1), Nfft(2), Nfft(3)];
total_size = prod(grid_size);

% 判断当前是否为“按需计算模式”（save_mem=true 时传入的 fft/idx 通常为空或非cell）
is_on_the_fly = isempty(fft) || ~iscell(fft);

%% Initialize FFT arrays
if use_gpu
    try
        if is_on_the_fly
            % 按需模式：直接使用单个 struct
            fft.Nfft1 = gpuArray.zeros(grid_size);
            fft.Nfft2 = gpuArray.zeros(grid_size);
        else
            % 预计算模式：存入 cell
            fft{iq}.Nfft1 = gpuArray.zeros(grid_size);
            fft{iq}.Nfft2 = gpuArray.zeros(grid_size);
        end
        use_gpu_actual = true;
    catch
        warning('GPU memory insufficient for FFT grids. Falling back to CPU.');
        use_gpu_actual = false;
        if is_on_the_fly
            fft.Nfft1 = zeros(grid_size);
            fft.Nfft2 = zeros(grid_size);
        else
            fft{iq}.Nfft1 = zeros(grid_size);
            fft{iq}.Nfft2 = zeros(grid_size);
        end
    end
else
    if is_on_the_fly
        fft.Nfft1 = zeros(grid_size);
        fft.Nfft2 = zeros(grid_size);
    else
        fft{iq}.Nfft1 = zeros(grid_size);
        fft{iq}.Nfft2 = zeros(grid_size);
    end
    use_gpu_actual = false;
end

if is_on_the_fly
    fft.size = total_size;
    fft.grid_size = grid_size;
    fft.use_gpu = use_gpu_actual;
else
    fft{iq}.size = total_size;
    fft{iq}.grid_size = grid_size;
    fft{iq}.use_gpu = use_gpu_actual;
end

%% gvec to fft grid
g_kq = wfnkq.mill;
g_k  = wfnk.mill;
g_q  = pol.mtx{iq};

bidx_kq = g2fft_grid(g_kq, Nfft(1), Nfft(2), Nfft(3));
bidx_k  = g2fft_grid(g_k,  Nfft(1), Nfft(2), Nfft(3));
bidx_q  = g2fft_grid(g_q,  Nfft(1), Nfft(2), Nfft(3));

%% Compute linear indices
lin_idx_kq = bidx_kq(:,1) + (bidx_kq(:,2)-1)*grid_size(1) + (bidx_kq(:,3)-1)*grid_size(1)*grid_size(2);
lin_idx_k  = bidx_k(:,1)  + (bidx_k(:,2)-1)*grid_size(1)  + (bidx_k(:,3)-1)*grid_size(1)*grid_size(2);
lin_idx_q  = bidx_q(:,1)  + (bidx_q(:,2)-1)*grid_size(1)  + (bidx_q(:,3)-1)*grid_size(1)*grid_size(2);

%% Store indices according to mode
if is_on_the_fly
    % === 按需计算模式（save_mem=true）：返回单个 struct ===
    idx.k  = lin_idx_k;
    idx.q  = lin_idx_q;
    idx.kq = lin_idx_kq;
    
    if use_gpu_actual
        idx.k  = gpuArray(idx.k);
        idx.q  = gpuArray(idx.q);
        idx.kq = gpuArray(idx.kq);
    end
else
    % === 预计算模式（save_mem=false）：存入 cell 数组 ===
    if use_gpu_actual
        idx.kq{iq, ik} = gpuArray(lin_idx_kq);
        idx.k{iq, ik}  = gpuArray(lin_idx_k);
        idx.q{iq}      = gpuArray(lin_idx_q);
    else
        idx.kq{iq, ik} = lin_idx_kq;
        idx.k{iq, ik}  = lin_idx_k;
        idx.q{iq}      = lin_idx_q;
    end
end

end