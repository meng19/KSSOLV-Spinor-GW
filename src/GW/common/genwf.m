function wfnkq = genwf_gpu(rkq, gr, gvec, syms, sys, options, cutoff, use_gpu)
% Check if GPU is requested and available
if use_gpu
    try
        % Test if GPU is available
        gpuArray();
        use_gpu_actual = true;
    catch
        warning('GPU requested but not available. Falling back to CPU.');
        use_gpu_actual = false;
    end
else
    use_gpu_actual = false;
end

% Precompute frequently used values
nspin = sys.nspin;
nspinor = sys.nspinor;
nbnd = sys.nbnd;

[ikrkq, itqq, kgqq] = find_kpt_match(gr, syms, rkq);
[ekin, isrtc_kq] = sortrx(rkq, gvec.ng, gvec.mill, sys);
nmtx = gcutoff(gvec.ng, ekin, isrtc_kq, cutoff);

kadd = g2fft_index(options.mill{1, ikrkq}, sys);
isortc = gvec.index_vec(kadd);

% Preallocate isorti for better performance
isorti = zeros(max(isortc(1:nmtx)), 1);
for i = 1:nmtx
    isorti(isortc(i)) = i;
end

ind = gmap(gvec, syms, nmtx, itqq, kgqq, isrtc_kq, isorti, sys);

%% assignment
wfnkq.ikq = ikrkq;
wfnkq.mill = gvec.mill(isrtc_kq(1:nmtx), :);
wfnkq.ind = ind;
wfnkq.psi = cell(nspin, nspinor);

% Pre-extract wavefunction data
wavefuncell_ikrkq = options.X0.wavefuncell(ikrkq, :);

% Process all spin and spinor components
for ispin = 1:nspin
    wavefunc_ispin = wavefuncell_ikrkq{ispin}.psi;
    
    for ispinor = 1:nspinor
        % Calculate indices for this spinor component
        start_idx = 1 + (ispinor - 1) * nmtx;
        end_idx = ispinor * nmtx;
        
        % Extract the wavefunction slice
        wfn_ispinor = wavefunc_ispin(start_idx:end_idx, :);
        
        % Move to GPU if requested and available
        if use_gpu_actual
            try
                wfn_ispinor = gpuArray(wfn_ispinor);
            catch ME
                if strcmp(ME.identifier, 'parallel:gpu:array:OOM')
                    warning('GPU memory insufficient. Falling back to CPU.');
                    use_gpu_actual = false;
                else
                    rethrow(ME);
                end
            end
        end
        
        % Apply indexing transformation
        wfnkq.psi{ispin, ispinor} = wfn_ispinor(ind, :);
    end
end

%% In spinor case, rotate spinors using optimized matrix operations
if nspinor == 2 && itqq ~= 1
    umtrx = susymmetries(sys.bvec, syms.mtrx{itqq}, itqq);
    
    % Move rotation matrix to GPU if using GPU
    if use_gpu_actual
        try
            umtrx = gpuArray(umtrx);
        catch
            % If GPU memory fails for umtrx, continue with CPU
        end
    end
    
    for ispin = 1:nspin
        % Extract both spinor components for this spin
        psi1 = wfnkq.psi{ispin, 1};
        psi2 = wfnkq.psi{ispin, 2};
        
        % OPTIMIZATION: Process ALL bands at once using matrix operations
        % Combine spinor components: [nmtx, nbnd, 2]
        spinor_combined = cat(3, psi1, psi2);
        
        if use_gpu_actual
            % GPU-accelerated batch processing using pagemtimes
            try
                % Reshape for efficient GPU computation: [nmtx, nbnd, 2] -> [nmtx, 2, nbnd]
                spinor_permuted = permute(spinor_combined, [1, 3, 2]);
                
                % Apply rotation to all bands simultaneously
                % umtrx is [2, 2], spinor_permuted is [nmtx, 2, nbnd]
                % Result: [nmtx, 2, nbnd] * [2, 2] -> [nmtx, 2, nbnd]
                rotated_spinor = pagemtimes(spinor_permuted, umtrx);
                
                % Reshape back: [nmtx, 2, nbnd] -> [nmtx, nbnd, 2]
                rotated_spinor_perm = permute(rotated_spinor, [1, 3, 2]);
                
            catch
                % Fall back to CPU if GPU batch operation fails
                use_gpu_actual = false;
                % Use CPU matrix multiplication for all bands
                rotated_spinor_perm = zeros(nmtx, nbnd, 2);
                for iband = 1:nbnd
                    rotated_spinor_perm(:, iband, :) = reshape(spinor_combined(:, iband, :) * umtrx, [nmtx, 1, 2]);
                end
            end
        else
            % CPU version: process all bands at once using matrix operations
            % Reshape to [nmtx * nbnd, 2] for batch multiplication
            spinor_reshaped = reshape(spinor_combined, [nmtx * nbnd, 2]);
            
            % Apply rotation to all bands
            rotated_reshaped = spinor_reshaped * umtrx;
            
            % Reshape back to original dimensions
            rotated_spinor_perm = reshape(rotated_reshaped, [nmtx, nbnd, 2]);
        end
        
        % Split back into spinor components
        psi1 = rotated_spinor_perm(:, :, 1);
        psi2 = rotated_spinor_perm(:, :, 2);
        
        % Update the cell array
        wfnkq.psi{ispin, 1} = psi1;
        wfnkq.psi{ispin, 2} = psi2;
        
        % Check norm (only on CPU to avoid GPU-CPU transfers)
        if ~use_gpu_actual
            wfn_ispin = [psi1; psi2];
            checknorm(wfn_ispin, ikrkq, ispin);
        end
    end
end

% Ensure all data is on CPU if GPU failed during processing
if use_gpu && ~use_gpu_actual
    for ispin = 1:nspin
        for ispinor = 1:nspinor
            if isa(wfnkq.psi{ispin, ispinor}, 'gpuArray')
                wfnkq.psi{ispin, ispinor} = gather(wfnkq.psi{ispin, ispinor});
            end
        end
    end
    
    % Re-check norms if we had to fall back to CPU
    if nspinor == 2 && itqq ~= 1
        for ispin = 1:nspin
            wfn_ispin = [wfnkq.psi{ispin, 1}; wfnkq.psi{ispin, 2}];
            checknorm(wfn_ispin, ikrkq, ispin);
        end
    end
end
end