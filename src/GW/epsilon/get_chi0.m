function chi0 = get_chi0(gme_q)
% gme_q形状为 [nmtx, nvbands, ncbands, nkpts]
[nmtx, nvbands, ncbands, nkpts] = size(gme_q);

% 重塑为矩阵 [nmtx, nvbands*ncbands*nkpts]
X = reshape(gme_q, nmtx, []);

% 在计算矩阵乘法: conj(X) * X.'
if ~isempty(X)
    chi0 = -conj(X) * X.';
end
end