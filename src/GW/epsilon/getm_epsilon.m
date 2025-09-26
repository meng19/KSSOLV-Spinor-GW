function gme = getm_epsilon_gpu(ivin, icin, wfnkq, wfnk, fft, idx, ispin, nspinor, options, use_gpu)
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

%% get eden and multiply
eval = options.ev(ivin, wfnkq.ikq, ispin);
econd = options.ev(icin, wfnk.ikq, ispin);
eden = 1/sqrt(eval-econd);

if use_gpu
    eden_gpu = gpuArray(eden);
    gme = gme * eden_gpu;
else
    gme = gme * eden;
end
end