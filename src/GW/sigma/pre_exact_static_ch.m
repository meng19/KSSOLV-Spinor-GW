function [igpp, valid_indices]= pre_exact_static_ch(fbz, gvec, indrk, iq, use_gpu)
% for an indrk(iq) kponit, calculate valid G_q(ig, :) - G_q(igp, :)
if use_gpu
    % GPU版本
    G_q = fbz.mtx_cutoff{1, indrk(iq)};
    ncouls = fbz.nmtx_cutoff(1, indrk(iq));
    
    % 将数据转移到GPU
    G_q_gpu = gpuArray(G_q);
    isorti_gpu = gpuArray(fbz.isorti(:, 1));
    
    % 创建网格索引
    [ig, igp] = meshgrid(1:ncouls, 1:ncouls);
    ig_gpu = gpuArray(ig(:));
    igp_gpu = gpuArray(igp(:));
    linearidx_gpu = sub2ind([ncouls ncouls], ig_gpu, igp_gpu);
    
    % 计算gpp = G_q(ig, :) - G_q(igp, :)
    G_q_ig = G_q_gpu(ig_gpu, :);
    G_q_igp = G_q_gpu(igp_gpu, :);
    gpp_gpu(linearidx_gpu, :) = G_q_ig - G_q_igp;
    
    gvec_gpu.mill = gpuArray(gvec.mill);
    gvec_gpu.index_vec = gpuArray(gvec.index_vec);
    gvec_gpu.nr = gvec.nr;
    gvec_gpu.ng = gvec.ng;
    gvec_gpu.n1 = gvec.n1;
    gvec_gpu.n2 = gvec.n2;
    gvec_gpu.n3 = gvec.n3;
    
    igpp_gpu = findvector(gpp_gpu, gvec_gpu);
    
    % 应用排序索引
    igpp = isorti_gpu(igpp_gpu);
    
    % 创建有效索引
    valid_indices = (igpp >= 1 & igpp <= fbz.nmtx(1, 1));
else
    G_q = fbz.mtx_cutoff{1, indrk(iq)};
    ncouls = fbz.nmtx_cutoff(1, indrk(iq));
    isorti = fbz.isorti(:, 1);
    
    [ig, igp] = meshgrid(1:ncouls, 1:ncouls);
    linearidx = sub2ind([ncouls ncouls], ig, igp);
    
    gpp = zeros(ncouls * ncouls, 3);
    gpp(linearidx, :) = G_q(ig, :) - G_q(igp, :);
    igpp = findvector(gpp, gvec);
    igpp = isorti(igpp);
    
    valid_indices = (igpp >= 1 & igpp <= fbz.nmtx(1, 1));
end
end