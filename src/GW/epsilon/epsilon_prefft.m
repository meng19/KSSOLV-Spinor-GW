function [fft, idx] = epsilon_prefft(wfnkq, wfnk, iq, ik, pol, fft, idx, use_gpu)
% For iq, precompute indices of wav in fft_grid in epsilon with GPU support
%% set up fft_grid
Nfac = 3;
Nfft = zeros(1, 3);
for i = 1:3
    Nfft(i) = pol.fftgrid{iq}(1, i);
    while ~check_FFT_size(Nfft(i), Nfac)
        Nfft(i) = Nfft(i) + 1;
    end
end

% Precompute grid size
grid_size = [Nfft(1), Nfft(2), Nfft(3)];
total_size = prod(grid_size);

% Initialize FFT arrays with appropriate data type and device
if use_gpu
    try
        fft{iq}.Nfft1 = gpuArray.zeros(grid_size, 'single');
        fft{iq}.Nfft2 = gpuArray.zeros(grid_size, 'single');
        use_gpu_actual = true;
    catch
        warning('GPU memory insufficient for FFT grids. Falling back to CPU.');
        fft{iq}.Nfft1 = zeros(grid_size, 'single');
        fft{iq}.Nfft2 = zeros(grid_size, 'single');
        use_gpu_actual = false;
    end
else
    fft{iq}.Nfft1 = zeros(grid_size, 'single');
    fft{iq}.Nfft2 = zeros(grid_size, 'single');
    use_gpu_actual = false;
end

fft{iq}.size = total_size;

%% gvec to fft grid, only translate
g_kq = wfnkq.mill;
g_k = wfnk.mill;
g_q = pol.mtx{iq};

% Compute FFT grid indices using vectorized operations
bidx_kq = g2fft_grid(g_kq, Nfft(1), Nfft(2), Nfft(3));
bidx_k = g2fft_grid(g_k, Nfft(1), Nfft(2), Nfft(3));
bidx_q = g2fft_grid(g_q, Nfft(1), Nfft(2), Nfft(3));

%% Precompute indices using optimized sub2ind
% Precompute linear indices for faster computation
if use_gpu_actual
    try
        % GPU version of sub2ind
        lin_idx_kq = bidx_kq(:,1) + (bidx_kq(:,2)-1)*grid_size(1) + (bidx_kq(:,3)-1)*grid_size(1)*grid_size(2);
        lin_idx_k = bidx_k(:,1) + (bidx_k(:,2)-1)*grid_size(1) + (bidx_k(:,3)-1)*grid_size(1)*grid_size(2);
        lin_idx_q = bidx_q(:,1) + (bidx_q(:,2)-1)*grid_size(1) + (bidx_q(:,3)-1)*grid_size(1)*grid_size(2);
        
        % Move indices to GPU
        idx.kq{iq, ik} = gpuArray(lin_idx_kq);
        idx.k{iq, ik} = gpuArray(lin_idx_k);
        idx.q{iq} = gpuArray(lin_idx_q);
    catch
        % Fall back to CPU if GPU memory is insufficient
        use_gpu_actual = false;
        idx.kq{iq, ik} = bidx_kq(:,1) + (bidx_kq(:,2)-1)*grid_size(1) + (bidx_kq(:,3)-1)*grid_size(1)*grid_size(2);
        idx.k{iq, ik} = bidx_k(:,1) + (bidx_k(:,2)-1)*grid_size(1) + (bidx_k(:,3)-1)*grid_size(1)*grid_size(2);
        idx.q{iq} = bidx_q(:,1) + (bidx_q(:,2)-1)*grid_size(1) + (bidx_q(:,3)-1)*grid_size(1)*grid_size(2);
    end
else
    % CPU version
    idx.kq{iq, ik} = bidx_kq(:,1) + (bidx_kq(:,2)-1)*grid_size(1) + (bidx_kq(:,3)-1)*grid_size(1)*grid_size(2);
    idx.k{iq, ik} = bidx_k(:,1) + (bidx_k(:,2)-1)*grid_size(1) + (bidx_k(:,3)-1)*grid_size(1)*grid_size(2);
    idx.q{iq} = bidx_q(:,1) + (bidx_q(:,2)-1)*grid_size(1) + (bidx_q(:,3)-1)*grid_size(1)*grid_size(2);
end

% Store grid dimensions for later use
fft{iq}.grid_size = grid_size;
fft{iq}.use_gpu = use_gpu_actual;
end