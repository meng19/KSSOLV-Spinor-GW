function eden = get_eden(ivin, icin, wfnkq, wfnk, ispin, options, freq, eps)
%% get eden with frequency dependence
eval = options.ev(ivin, wfnkq.ikq, ispin);
econd = options.ev(icin, wfnk.ikq, ispin);

if eps.freq_dep == 2
    % 频率依赖的介电函数计算
    energy_diff = eval - econd;
    
    % 计算频率依赖的情况
    % denominator = 2 / (1 / (freq + energy_diff) + 1 / (-freq + energy_diff));
    % denominator = (energy_diff^2 - freq^2) / energy_diff;
    % eden = 1 / sqrt(denominator); 使用sqrt乘以每个gme会造成误差？
    eden = energy_diff / (energy_diff^2 - freq^2);
elseif eps.freq_dep == 0
    eden = 1/(eval-econd);
end