function gme = isdf_sigma_matrix_elements_from_real(state_real, sum_real, idx_q, fftgrid, options)
%ISDF_SIGMA_MATRIX_ELEMENTS_FROM_REAL Approximate sigma matrix elements.
%   Computes FFT_G(conj(state_real) .* sum_real(:, nn)) for all nn.

phi = conj(state_real(:));
gme3 = isdf_epsilon_matrix_elements_from_real(phi, sum_real, idx_q, fftgrid, options);
gme = reshape(gme3, length(idx_q), size(sum_real, 2));
end
