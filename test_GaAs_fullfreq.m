% Cleaning up the workspace and command line window:
clc
clear all;
close all;
randn('state', 0);
rand('state', 0);

% 启动 Profiler
profile on

% Initializing the environment path for KSSOLV:
KSSOLV_startup;

read_vxc = 1;
[sys, options, syms] = read_qe_gw('./example/qe_data', read_vxc);
[sys, options] = gw_setup(sys, options);

eps.nbnd = 999;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
eps.freq_dep = 2;
eps.freq_dep_method = 2;
eps.freq_cutoff = 200;
eps.delta_freq = 15;
eps.nfreq_imag = 15;
omega = 0;
eta = 0;
eps.cutoff = 20; % Ry
eps.coul_cut = 'spherical_truncation';
eps.coul_cutoff = 5; %coulomb truncation radius in epsilon
eps.use_gpu = 1;
eps.save_mem = 0;
eps.precompute_wav = 0;
eps = epsilon(sys, options, syms, eps);
save('epsilon_spin.mat','options','syms','sys','eps','-v7.3')

sig.nbnd = 999;
sig.ndiag_min = 37;
sig.ndiag_max = 52;
sig.freq_dep = 2;
sig.freq_dep_method = 2;
sig.freq_grid_shift = 2;
sig.max_freq_eval = 2;
sig.delta_freq_eval = 0.2;
sig.cd_int_method = 0;
sig.coul_cut = 'spherical_truncation';
sig.coul_cutoff = 5; %coulomb truncation radius in sigma
sig.no_symmetries_q_grid = 0;
sig.exact_static_ch = 1;
sig.use_gpu = 1;
sig.precompute_wav = 0;
sig = sigma(eps, sig, sys, options, syms);

fid = fopen('eqp0.txt', 'w');
fprintf(fid, [repmat('%.15g\t', 1, size(sig.eqp0, 2)-1), '%.15g\n'], sig.eqp0.');
fclose(fid);

fid = fopen('eqp1.txt', 'w');
fprintf(fid, [repmat('%.15g\t', 1, size(sig.eqp1, 2)-1), '%.15g\n'], sig.eqp1.');
fclose(fid);
% 停止 Profiler 并生成报告
profile off
% 保存 Profiler 数据
p = profile('info');
profsave(p, 'profile');

save('sigma_spin.mat','options','syms','sys','sig','-v7.3')
