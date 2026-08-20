function screened = sigma_gpu_screened_w(screened)
%SIGMA_GPU_SCREENED_W Move screened-W numeric fields to GPU.

if isempty(screened)
    return;
end
fields = fieldnames(screened);
for ifield = 1:numel(fields)
    value = screened.(fields{ifield});
    if isnumeric(value) && ~isa(value, 'gpuArray')
        screened.(fields{ifield}) = gpuArray(value);
    end
end
end
