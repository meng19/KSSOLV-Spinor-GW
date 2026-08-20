function values = randn_like(nrows, ncols, reference)
%RANDN_LIKE Create random normal values matching the reference precision/device.

if isa(reference, 'gpuArray')
    values = randn(nrows, ncols, 'like', real(reference));
else
    values = randn(nrows, ncols);
    if isa(reference, 'single')
        values = single(values);
    end
end
end
