function values = real_component( ...
    wfn, fft_template, idx, ispin, ispinor, band_list, progress)
%ISDF.REAL_COMPONENT Convert one spinor component to real space.

if nargin < 7
    progress = [];
end
ngrid = numel(fft_template);
values = complex(zeros(ngrid, numel(band_list), 'like', fft_template));
for iband = 1:numel(band_list)
    fft_box = fft_template;
    fft_box(idx) = wfn.psi{ispin, ispinor}(:, band_list(iband));
    band_real = fftn(fft_box);
    values(:, iband) = band_real(:);
    if isa(progress, 'function_handle')
        progress(iband, numel(band_list));
    end
end
end
