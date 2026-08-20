function kernel = page_project_kernel(left_projector, k_mu, right_projector)
%PAGE_PROJECT_KERNEL Project frequency-page reduced screened kernels.

nfreq = size(k_mu, 3);

if exist('pagemtimes', 'file') == 2 || exist('pagemtimes', 'builtin') == 5
    try
        kernel = pagemtimes(pagemtimes(left_projector, k_mu), ...
            right_projector);
        return;
    catch ME
        if ~any(strcmp(ME.identifier, { ...
                'MATLAB:pagemtimes:innerdim', ...
                'MATLAB:pagemtimes:InputSizeMismatch'}))
            rethrow(ME);
        end
    end
end

kernel = local_project_with_gemm(left_projector, k_mu, right_projector);
end

function kernel = local_project_with_gemm(left_projector, k_mu, right_projector)
% Use two larger GEMMs instead of nfreq small GEMMs.

nleft = size(left_projector, 1);
nmu = size(k_mu, 1);
nfreq = size(k_mu, 3);
nright = size(right_projector, 2);

left_k = left_projector * reshape(k_mu, nmu, nmu * nfreq);
left_k_pages = reshape(left_k, nleft, nmu, nfreq);
left_k_rows = reshape(permute(left_k_pages, [1, 3, 2]), ...
    nleft * nfreq, nmu);
kernel_rows = left_k_rows * right_projector;
kernel = permute(reshape(kernel_rows, nleft, nfreq, nright), ...
    [1, 3, 2]);
end
