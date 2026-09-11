function value = gw_gather_if_gpu(value)
%GW_GATHER_IF_GPU Gather gpuArray values and otherwise return them unchanged.

if isa(value, 'gpuArray')
    value = gather(value);
end
end
