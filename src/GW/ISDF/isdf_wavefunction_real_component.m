function real_values = isdf_wavefunction_real_component(wfn, fft_template, idx, ispin, ispinor, band_list)
%ISDF_WAVEFUNCTION_REAL_COMPONENT Convert one spinor component to real space.

ngrid = numel(fft_template);
real_values = zeros(ngrid, length(band_list));

for ib = 1:length(band_list)
    fft_box = fft_template;
    fft_box(idx) = wfn.psi{ispin, ispinor}(:, band_list(ib));
    band_real = fftn(fft_box);
    real_values(:, ib) = band_real(:);
end
end
