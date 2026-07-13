function gme_batch = isdf_sigma_batch(wfnkq, wfnk, fft, idx, ispin, nspinor, state_band, sum_bands, isdf_options)
%ISDF_SIGMA_BATCH Compute sigma matrix elements for one state and many bands.

nq = length(idx.q);
nbands = length(sum_bands);
gme_batch = zeros(nq, nbands);

fftgrid = size(fft.Nfft1);
for ispinor = 1:nspinor
    state_real = isdf_wavefunction_real_component( ...
        wfnk, fft.Nfft1, idx.k, ispin, ispinor, state_band);
    sum_real = isdf_wavefunction_real_component( ...
        wfnkq, fft.Nfft2, idx.kq, ispin, ispinor, sum_bands);

    gme_batch = gme_batch + isdf_sigma_matrix_elements_from_real( ...
        state_real, sum_real, idx.q, fftgrid, isdf_options);
end
end
