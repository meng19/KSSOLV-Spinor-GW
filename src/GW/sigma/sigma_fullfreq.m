function [asx_loc, ach_loc, achx_loc, omega, iw_lda] = sigma_fullfreq(asx_loc, ach_loc, in, nn, ik, iq, occ_kq, emf, ispin, aqsn, aqsm, eps_inv_I_coul, sig)
% sigma_ch - contour-deformation 相关自能计算
% 输入：
%   in, nn          带索引（通常 in 是 n1，nn 是 n1true 或目标带）
%   emf             准粒子能量 [nbands × nspin]
%   ispin           自旋索引
%   sig             参数结构体（包含频率网格、方法等）
%   aqsntemp        <ψ_nk | e^{iGr} | ψ_m> 矩阵 (ncouls × n1max)
%   aqsmtemp        对应的共轭或另一形式 (ncouls × n1max)


ryd = 13.6056923;  % Rydberg 常数 (eV)

e_nk   = emf(in, ik, ispin) * ryd;    % ε_{n k}
e_n1kq = emf(nn, iq, ispin) * ryd;    % ε_{n' k+q}
aqsmn = aqsm .* aqsn';
M_flat = aqsmn(:);

% 确定实频率网格起点
if sig.freq_grid_shift < 2
    freq0 = sig.min_freq_eval;
else
    freq0 = e_nk - sig.delta_freq_eval * (sig.nfreq_grid - 1) / 2;
end

% 实频率网格点 ω_i
nfreqeval = sig.nfreq_grid;
iw        = 0 : nfreqeval-1;
omega     = freq0 + iw * sig.delta_freq_eval;
wxi       = omega - e_n1kq;           % ω_i - ε_{n'k+q}
wxi = wxi';
wxi(wxi == 0) = 1e-6;

% 找到最接近 e_nk 的网格点（常用于对齐 LDA 能量）
[~, iw_lda] = min(abs(omega - e_nk));

% 常用参数提取
nfreq_real = sig.nfreq_integral_real;
nfreq_imag = sig.nfreq_integral_imag;
ncouls     = size(aqsn, 1);

% 初始化
sres = zeros(nfreqeval, 1) + 0i;   % residue 贡献
sint = zeros(nfreqeval, 1) + 0i;   % 虚频率积分贡献

% ───────────────────────────────────────────────
% 1. 静态 COH 贡献（如果开启 exact_ch）
% ───────────────────────────────────────────────
if isfield(sig, 'exact_static_ch') && sig.exact_static_ch == 1
%     achx_loc = -0.5 * sum(sum(aqsm .* aqsn' .* eps_inv_I_coul(:,:,1) * 0.5));
    E_1 = reshape(eps_inv_I_coul(:,:,1), ncouls*ncouls, 1);
    achx_loc = -0.25 * M_flat.' * E_1;
else
    achx_loc = 0;
end

% ───────────────────────────────────────────────
% 2. Residue 部分（实频率轴上的极点贡献）
% ───────────────────────────────────────────────
% 判断哪些 ω_i 需要 residue 贡献
need_res = ((wxi >= 0) ~= (occ_kq > 0));
idx_need = find(need_res);

if ~isempty(idx_need)
    wx_abs_need = abs(wxi(idx_need));
    if nfreq_real <= 1
        sres_omega_need = sum(aqsmn .* eps_inv_I_coul(:,:,1), 'all') * ones(length(idx_need), 1);
    else
        ifreq_need = zeros(length(idx_need), 1, 'int32');
        for k = 1:nfreq_real-1
            mask = wx_abs_need >= sig.freq_integral(k) & wx_abs_need < sig.freq_integral(k+1);
            ifreq_need(mask) = k;
        end
        ifreq_need(ifreq_need == 0) = nfreq_real - 1;
        
        f_low  = sig.freq_integral(ifreq_need)';
        f_high = sig.freq_integral(ifreq_need + 1)';
        df     = f_high - f_low;
        fact2  = (wx_abs_need - f_low) ./ df;
        fact1  = 1 - fact2;
        
        sres_omega_need = zeros(length(idx_need), 1);
        
        for i = 1:length(idx_need)
            k = ifreq_need(i);
            E_k  = eps_inv_I_coul(:,:,k);     % 列向量
            E_kp = eps_inv_I_coul(:,:,k+1);
            
            % 向量内积
            sres_omega_need(i) = fact1(i) * (M_flat.' * E_k(:)) + ...
                fact2(i) * (M_flat.' * E_kp(:));
        end
    end
    
    % 填回 sres（只更新需要的点）
    occ_sign_vec = sign(wxi(idx_need));
    occ_sign_vec(occ_sign_vec == 0) = 1;
    sres(idx_need) = occ_sign_vec .* sres_omega_need;
end
% ───────────────────────────────────────────────
% 3. 虚频率积分部分（最常用 method=0）
% ───────────────────────────────────────────────
% sW_imag = sum(sum(aqsm .* aqsn' .* eps_inv_I_coul(:,:,nfreq_real+1 : nfreq_real+nfreq_imag)));
E_flat = reshape(eps_inv_I_coul(:,:, nfreq_real+1 : nfreq_real + nfreq_imag), ncouls*ncouls, nfreq_imag);
sW_imag = M_flat.' * E_flat;
sW_imag = -sW_imag(:);

imag_freqs = imag(sig.freq_integral(nfreq_real+1 : nfreq_real+nfreq_imag));

if ~isfield(sig, 'cd_int_method') || sig.cd_int_method == 0
    % method 0：简单中点梯形 + atan 积分
    sint = zeros(nfreqeval, 1) + 0i;
    
    % 第一段
    sW = sW_imag(1);
    fStart = imag_freqs(1);
    fMid   = (imag_freqs(1) + imag_freqs(2)) / 2;
    sint = sint + sW * (atan(fMid ./ wxi) - atan(fStart ./ wxi));
    
    % 中间段
    for iq = 2 : nfreq_imag-1
        sW = sW_imag(iq);
        fStart = (imag_freqs(iq-1) + imag_freqs(iq)) / 2;
        fMid   = (imag_freqs(iq)   + imag_freqs(iq+1)) / 2;
        sint = sint + sW * (atan(fMid ./ wxi) - atan(fStart ./ wxi));
    end
    
    % 最后一段（外推）
    sW = sW_imag(end);
    fStart = (imag_freqs(end-1) + imag_freqs(end)) / 2;
    fEnd   = (-imag_freqs(end-1) + 3*imag_freqs(end)) / 2;
    sint = sint + sW * (atan(fEnd ./ wxi) - atan(fStart ./ wxi));
else
    warning('本版本只实现了 cd_int_method = 0，其他积分方法待补充');
end

% ───────────────────────────────────────────────
% 最终自能
% ───────────────────────────────────────────────
asx_loc = asx_loc + sres;
ach_loc = ach_loc + sint / pi;

% 如果需要单独保存 LDA 点处的 achtD 值，可以这样：
% achtD_n1(nn) = ach(iw_lda);   % 需要在外部定义 achtD_n1 数组

end