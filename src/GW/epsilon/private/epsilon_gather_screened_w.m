function screened = epsilon_gather_screened_w(screened)
%EPSILON_GATHER_SCREENED_W Gather gpuArray fields in screened-W structs.

fields = fieldnames(screened);
for ifield = 1:numel(fields)
    value = screened.(fields{ifield});
    if isa(value, 'gpuArray')
        screened.(fields{ifield}) = gather(value);
    end
end
end
