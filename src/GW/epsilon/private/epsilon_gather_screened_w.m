function screened = epsilon_gather_screened_w(screened)
%EPSILON_GATHER_SCREENED_W Gather gpuArray fields in screened-W structs.
%   On CPU, return the original struct without structfun rebuilding a
%   second output struct while the large reduced-W pages are live.

fields = fieldnames(screened);
has_gpu_field = false;
for ifield = 1:numel(fields)
    if isa(screened.(fields{ifield}), 'gpuArray')
        has_gpu_field = true;
        break;
    end
end
if ~has_gpu_field
    return;
end
for ifield = 1:numel(fields)
    field = fields{ifield};
    screened.(field) = gw_gather_if_gpu(screened.(field));
end
end
