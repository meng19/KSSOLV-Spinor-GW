function screened = epsilon_gather_screened_w(screened)
%EPSILON_GATHER_SCREENED_W Gather gpuArray fields in screened-W structs.

screened = structfun(@gw_gather_if_gpu, screened, ...
    'UniformOutput', false);
end
