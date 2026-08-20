function value = gather_if_gpu(value)
%GATHER_IF_GPU Return CPU data when the input is a gpuArray.

if isa(value, 'gpuArray')
    value = gather(value);
end
end
