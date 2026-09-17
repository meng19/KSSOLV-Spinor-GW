function ind_mu = qrcp_sample(products, rank_mu)
%QRCP_SAMPLE Select interpolation points from product rows.

try
    [~, ~, pivots] = qr(products', 0);
    
    if (0)
    % 用 SVD 右奇异向量做列选择（近似 CPQR 的前 n 个主元）
    [~, ~, V] = svds(products', rank_mu);
    % 杠杆分数：每列在主子空间中的"重要性"
    lev = sum(V.^2, 2);
    [~, pivots] = maxk(lev, rank_mu);   % 近似的前 n 个主元列
    end
catch ME
    if ~isa(products, 'gpuArray')
        rethrow(ME);
    end
    [~, ~, pivots] = qr(gather_if_gpu(products'), 0);
end
ind_mu = gather_if_gpu(reshape(pivots(1:rank_mu), 1, []));
end
