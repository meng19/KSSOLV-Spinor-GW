function [asx_loc, ach_loc, iw_lda] = sigma_fullfreq(asx_loc, ach_loc, in, nn, occ_kq, emf, ispin, aqsn, aqsm, eps_inv_I_coul, sig)
% sigma_ch - contour-deformation 相关自能计算
% 输入：
%   in, nn          带索引（通常 in 是 n1，nn 是 n1true 或目标带）
%   emf             准粒子能量 [nbands × nspin]
%   ispin           自旋索引
%   sig             参数结构体（包含频率网格、方法等）
%   aqsntemp        <ψ_nk | e^{iGr} | ψ_m> 矩阵 (ncouls × n1max)
%   aqsmtemp        对应的共轭或另一形式 (ncouls × n1max)


    ryd = 13.6056923;  % Rydberg 常数 (eV)

    e_nk   = emf(in, ispin) * ryd;    % ε_{n k}
    e_n1kq = emf(nn, ispin) * ryd;    % ε_{n' k+q}

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
    wxi(wxi == 0) = 1e-6; % 避免被除

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
        achx_loc = sum(sum(aqsm .* aqsn' .* eps_inv_I_coul(:,:,1) * 0.5));
    end

    % ───────────────────────────────────────────────
    % 2. Residue 部分（实频率轴上的极点贡献）
    % ───────────────────────────────────────────────
    for iw = 1:nfreqeval
        wx = wxi(iw);
        if (wx >=0) == (occ_kq > 0)
            continue
        end
        if wx >= 0
            occ_sign = 1;
        else
            occ_sign = -1;
        end
        wx_abs = abs(wx);

        % 找到 wx 所在的频率区间
        ifreq = 0;
        for k = 1:nfreq_real-1
            if wx_abs >= sig.freq_integral(k) && wx_abs < sig.freq_integral(k+1)
                ifreq = k;
                break;
            end
        end
        if ifreq == 0
            ifreq = nfreq_real - 1;  % 超出范围时取最后一个区间
        end

        % 线性插值 ε^{-1}(ω)
        if nfreq_real > 1
            df    = sig.freq_integral(ifreq+1) - sig.freq_integral(ifreq);
            fact1 = (sig.freq_integral(ifreq+1) - wx_abs) / df;
            fact2 = (wx_abs - sig.freq_integral(ifreq))   / df;
            eps_inv_I_coul_freq = fact1 * eps_inv_I_coul(:,:,ifreq) + fact2 * eps_inv_I_coul(:,:,ifreq+1);
        else
            eps_inv_I_coul_freq = eps_inv_I_coul(:,:,1);
        end

        % 计算 Σ_res ~ Σ_G G' <n| e^{iGr} |n'> ε^{-1}_{GG'}(ω) <n'| e^{-iG'r} |n> v_G
        sres_omega = sum(sum(aqsm .* aqsn' .* eps_inv_I_coul_freq));

        sres(iw) = occ_sign * sres_omega;
    end

    % ───────────────────────────────────────────────
    % 3. 虚频率积分部分（最常用 method=0）
    % ───────────────────────────────────────────────
    sW_imag = sum(sum(aqsm .* aqsn' .* eps_inv_I_coul(:,:,nfreq_real+1 : nfreq_real+nfreq_imag)));
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