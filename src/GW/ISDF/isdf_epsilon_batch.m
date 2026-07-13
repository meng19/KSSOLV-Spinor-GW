function gme_batch = isdf_epsilon_batch(wfnkq, wfnk, fft, idx, ispin, nspinor, valence_bands, conduction_bands, isdf_options)
%ISDF_EPSILON_BATCH Compute a batch of epsilon matrix elements with ISDF.

nq = length(idx.q);
nv = length(valence_bands);
nc = length(conduction_bands);
gme_batch = zeros(nq, nv, nc);

fftgrid = size(fft.Nfft1);
for ispinor = 1:nspinor
    valence_real = isdf_wavefunction_real_component( ...
        wfnkq, fft.Nfft1, idx.kq, ispin, ispinor, valence_bands);
    conduction_real = isdf_wavefunction_real_component( ...
        wfnk, fft.Nfft2, idx.k, ispin, ispinor, conduction_bands);

    phi = conj(valence_real);
    psi = conduction_real;
    gme_batch = gme_batch + isdf_epsilon_matrix_elements_from_real( ...
        phi, psi, idx.q, fftgrid, isdf_options);
end
end
