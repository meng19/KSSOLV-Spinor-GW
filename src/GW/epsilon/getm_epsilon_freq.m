function [gme, eden] = getm_epsilon_freq(ivin, icin, wfnkq, wfnk, fft, idx, ispin, nspinor, use_gpu)
% 频率依赖的矩阵元计算函数
% 输入参数：
%   ivin: 价带索引
%   icin: 导带索引
%   wfnkq: k+q点的波函数
%   wfnk: k点的波函数
%   fft: FFT网格信息
%   idx: 索引信息
%   ispin: 自旋索引
%   nspinor: 自旋轨道数
%   options: 选项参数
%   freq: 频率值
%   use_gpu: 是否使用GPU

if use_gpu
    % GPU FFT优化
    persistent gpu_fft_initialized
    if isempty(gpu_fft_initialized)
        % GPU FFT预热
        temp = gpuArray(complex(zeros(32, 32, 32, 'single')));
        fftn(temp);
        gpu_fft_initialized = true;
    end
else
    % FFTW optimization for CPU
    persistent fftw_planned
    if isempty(fftw_planned)
        fftw('planner', 'measure');
        fftw_planned = true;
    end
end

% Preallocate grid
Nfft1 = fft.Nfft1;
Nfft2 = fft.Nfft2;

%% Main computation loop
for ispinor = 1:nspinor
    % 执行填充操作
    Nfft1(idx.kq) = wfnkq.psi{ispin, ispinor}(:, ivin);
    Nfft2(idx.k) = wfnk.psi{ispin, ispinor}(:, icin);
    
    % FFT计算
    fft_Nfft1 = fftn(Nfft1);
    fft_Nfft2 = fftn(Nfft2);
    
    % 计算卷积：conj(fft_Nfft1) .* fft_Nfft2
    conv = conj(fft_Nfft1) .* fft_Nfft2;
    
    % FFT并归一化
    result = fftn(conv) / fft.size;
    
    % 提取结果
    gme(:, ispinor) = result(idx.q);
end

%% Sum over spinor components
gme = sum(gme, 2);
end