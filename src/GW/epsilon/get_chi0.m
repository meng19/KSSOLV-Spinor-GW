function chi0 = get_chi0(gme_q, eden_q)
% get_chi0 计算 chi0
% 输入:
%   gme_q     [nmtx, nvbands, ncbands, nkpts]
%   eden_q    [nvbands, ncbands, nkpts, nfreq] （可选，缺省为 1）
% 输出:
%   chi0      [nmtx, nmtx, nfreq]

[nmtx, nvbands, ncbands, nkpts] = size(gme_q);
[nvbands, ncbands, nkpts, nfreq] = size(eden_q);

% 如果没有传入 eden_q 或为空，则设为 1
if nargin < 2 || isempty(eden_q)
    eden_q = 1;
end

% 重塑 gme_q 为 [nmtx, nvbands*ncbands*nkpts]
X = reshape(gme_q, nmtx, []);

% 处理 eden_q
if isscalar(eden_q) && eden_q == 1
    % 最常见情况：eden=1，直接计算 conj(X) * X.'
    if ~isempty(X)
        chi0 = -conj(X) * X.'; % 负号来自于共轭，兼容老版本，此时gme被乘以一个sqrt(-1)了?
    else
        chi0 = zeros(nmtx, nmtx);
    end
    
else
    % eden_q 不是标量1的情况，按原逻辑处理
    % 重塑为行向量 [1, nvbands*ncbands*nkpts]
    eden = reshape(eden_q, [], nfreq);
    
    if ~isempty(X)
        for ifreq = 1:nfreq
            chi0(:,:, ifreq) = conj(X) .* eden(:, ifreq)' * X.';
        end
    else
        chi0 = zeros(nmtx, nmtx, nfreq);
    end
end

end