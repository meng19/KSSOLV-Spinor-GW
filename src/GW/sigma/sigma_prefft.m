function idx = sigma_prefft(wfnkq, wfnk, q_mtx, iq, ik, sys, idx, use_gpu)
% For iq, precompute indices of wav in fft_grid in sigma with GPU support
%% set up fft_grid
Nfft = [sys.n1, sys.n2, sys.n3];
grid_size = [Nfft(1), Nfft(2), Nfft(3)];

%% gvec to fft grid, only translate
g_kq = wfnkq.mill;
g_k = wfnk.mill;
g_q = q_mtx;

% Vectorized grid index computation
bidx_kq = g2fft_grid(g_kq, Nfft(1), Nfft(2), Nfft(3));
bidx_k = g2fft_grid(g_k, Nfft(1), Nfft(2), Nfft(3));
bidx_q = g2fft_grid(g_q, Nfft(1), Nfft(2), Nfft(3));

%% Precompute indices using optimized linear index calculation
% Manual linear index calculation (faster than sub2ind)
if use_gpu
    try
        % GPU version
        lin_idx_kq = bidx_kq(:,1) + (bidx_kq(:,2)-1)*grid_size(1) + (bidx_kq(:,3)-1)*grid_size(1)*grid_size(2);
        lin_idx_k = bidx_k(:,1) + (bidx_k(:,2)-1)*grid_size(1) + (bidx_k(:,3)-1)*grid_size(1)*grid_size(2);
        lin_idx_q = bidx_q(:,1) + (bidx_q(:,2)-1)*grid_size(1) + (bidx_q(:,3)-1)*grid_size(1)*grid_size(2);
        
        % Move to GPU
        idx.kq{iq, ik} = gpuArray(lin_idx_kq);
        idx.k{ik} = gpuArray(lin_idx_k);
        idx.q{iq, ik} = gpuArray(lin_idx_q);
        
    catch ME
        if strcmp(ME.identifier, 'parallel:gpu:array:OOM')
            warning('GPU memory insufficient for sigma indices. Falling back to CPU.');
            % CPU fallback
            idx.kq{iq, ik} = bidx_kq(:,1) + (bidx_kq(:,2)-1)*grid_size(1) + (bidx_kq(:,3)-1)*grid_size(1)*grid_size(2);
            idx.k{ik} = bidx_k(:,1) + (bidx_k(:,2)-1)*grid_size(1) + (bidx_k(:,3)-1)*grid_size(1)*grid_size(2);
            idx.q{iq, ik} = bidx_q(:,1) + (bidx_q(:,2)-1)*grid_size(1) + (bidx_q(:,3)-1)*grid_size(1)*grid_size(2);
        else
            rethrow(ME);
        end
    end
else
    % CPU version
    idx.kq{iq, ik} = bidx_kq(:,1) + (bidx_kq(:,2)-1)*grid_size(1) + (bidx_kq(:,3)-1)*grid_size(1)*grid_size(2);
    idx.k{ik} = bidx_k(:,1) + (bidx_k(:,2)-1)*grid_size(1) + (bidx_k(:,3)-1)*grid_size(1)*grid_size(2);
    idx.q{iq, ik} = bidx_q(:,1) + (bidx_q(:,2)-1)*grid_size(1) + (bidx_q(:,3)-1)*grid_size(1)*grid_size(2);
end

% Store additional information for potential reuse
if ~isfield(idx, 'grid_size')
    idx.grid_size = grid_size;
end
if ~isfield(idx, 'use_gpu')
    idx.use_gpu = use_gpu && isa(idx.kq{iq, ik}, 'gpuArray');
end
end