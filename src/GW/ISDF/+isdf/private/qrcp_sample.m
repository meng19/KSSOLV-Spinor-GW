function ind_mu = qrcp_sample(products, rank_mu)
%QRCP_SAMPLE Select interpolation points from product rows.

try
    [~, ~, pivots] = qr(products', 0);
catch ME
    if ~isa(products, 'gpuArray')
        rethrow(ME);
    end
    [~, ~, pivots] = qr(gather_if_gpu(products'), 0);
end
ind_mu = gather_if_gpu(reshape(pivots(1:rank_mu), 1, []));
end
