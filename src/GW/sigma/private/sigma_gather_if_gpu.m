function value = sigma_gather_if_gpu(value)
%SIGMA_GATHER_IF_GPU Gather gpuArray values before storing sigma outputs.

if isa(value, 'gpuArray')
    value = gather(value);
end
end
