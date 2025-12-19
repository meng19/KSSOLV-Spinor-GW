function ws_dist = calculate_ws_distance(grid_size, metric_tensor, trunc_shift, n_cell)
% CALCULATE_WS_DISTANCE 计算在度规 metric_tensor 下，
% 网格点到其周期性镜像点的最小 Wigner-Seitz 距离。
% 参数说明：
% grid_size: [N1, N2, N3]
% metric_tensor: 3x3 矩阵 (adot)
% trunc_shift: [s1, s2, s3] 偏移量
% n_cell: 周期性扩展层数

N1 = grid_size(1); N2 = grid_size(2); N3 = grid_size(3);

% 1. 生成基础网格坐标 (基于 N1, N2, N3 顺序)
% grid_X, grid_Y, grid_Z 的尺寸均为 [N1, N2, N3]
[grid_X, grid_Y, grid_Z] = ndgrid(...
    (0:N1-1) + trunc_shift(1), ...
    (0:N2-1) + trunc_shift(2), ...
    (0:N3-1) + trunc_shift(3));

% 将网格拉直为 3 x (N1*N2*N3)
% 此时拉直的顺序是先 N1，再 N2，最后 N3
P_base = [grid_X(:), grid_Y(:), grid_Z(:)]';

% 2. 生成周期性位移向量 (L_scaled)
l_range = (-n_cell + 1) : n_cell;
[L1, L2, L3] = ndgrid(l_range, l_range, l_range);
shift_offsets = [L1(:)' * N1; L2(:)' * N2; L3(:)' * N3];

% 3. 初始化最小距离
num_points = size(P_base, 2);
min_dist_vec = inf(1, num_points, 'like', metric_tensor);

% 4. 向量化循环：遍历位移向量
% 对于每一个可能的周期性平移，计算所有网格点到该平移后的距离
for idx = 1:size(shift_offsets, 2)
    % 相对位置 X = P - L
    X = P_base - shift_offsets(:, idx);
    
    % 二次型计算距离平方 d^2 = X' * A * X
    dist_sq = sum(X .* (metric_tensor * X), 1);
    
    % 计算距离并更新最小值
    % 使用 max(..., 0) 确保数值稳定性（防止微小负数导致 sqrt 产生复数）
    min_dist_vec = min(min_dist_vec, sqrt(max(dist_sq, 0)));
end

% 5. 还原形状为 [N1, N2, N3]
ws_dist = reshape(min_dist_vec, [N1, N2, N3]);
end