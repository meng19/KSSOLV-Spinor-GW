% Validate optional ISDF epsilon path for MoS2 against the existing direct path.
clc;
clear all;
close all;
randn('state', 0);
rand('state', 0);

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(script_dir)));
addpath(repo_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(repo_root);
KSSOLV_startup;

read_vxc = 1;
[sys, options, syms] = read_qe_gw('.\example\qe_data\mos2_222_spinor', read_vxc);
[sys, options] = gw_setup(sys, options);

eps_base.nbnd = 29;
eps_base.nv = options.nv;
eps_base.nc = eps_base.nbnd - eps_base.nv;
eps_base.freq_dep = 0;
eps_base.freq_dep_method = 2;
eps_base.freq_cutoff = 200;
eps_base.delta_freq = 15;
eps_base.nfreq_imag = 15;
eps_base.cutoff = 10;
eps_base.coul_cut = 'spherical_truncation';
eps_base.coul_cutoff = 10;
eps_base.use_gpu = 0;
eps_base.save_mem = 0;
eps_base.precompute_wav = 0;

tic
eps_direct = epsilon(sys, options, syms, eps_base);
toc

eps_isdf_input = eps_base;
eps_isdf_input.isdf.enable = true;
eps_isdf_input.isdf.algorithm = 'matrix_elements';
eps_isdf_input.isdf.sample_method = 'qrcp';
eps_isdf_input.isdf.rank = 78;
eps_isdf_input.isdf.seed = 0;

tic
eps_isdf = epsilon(sys, options, syms, eps_isdf_input);
toc

max_relative_error = 0;
for iq = 1:length(eps_direct.inv)
    numerator = norm(eps_direct.inv{iq}(:) - eps_isdf.inv{iq}(:));
    denominator = max(1, norm(eps_direct.inv{iq}(:)));
    max_relative_error = max(max_relative_error, numerator / denominator);
end

fprintf('MoS2 epsilon ISDF validation max relative error = %.3e\n', max_relative_error);
assert(max_relative_error < 1e-8, ...
    'MoS2 epsilon ISDF validation failed: relative error %.3e', max_relative_error);
