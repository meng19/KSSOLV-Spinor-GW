function [eps, sig, sys, options, syms] = gw_run_json(config_file)
%GW_RUN_JSON Run a GW calculation described by a portable JSON file.
%   GW_RUN_JSON(CONFIG_FILE) may be called from any working directory.  All
%   relative paths in CONFIG_FILE are interpreted relative to that file.
%
%   Required JSON fields are:
%     qe_path, epsilon, sigma
%   Optional fields are:
%     read_vxc (true), rng_seed, save_file
%
%   Example launcher (the batch script can live anywhere):
%     matlab -batch "addpath('C:/path/KSSOLV-Spinor-GW'); \
%         gw_run_json('C:/path/case.json')"

if nargin ~= 1 || ~(ischar(config_file) || isstring(config_file))
    error('gw_run_json:InvalidInput', 'Provide exactly one JSON configuration file.');
end

config_file = char(config_file);
if ~isfile(config_file)
    error('gw_run_json:ConfigNotFound', 'Configuration file not found: %s', config_file);
end
config_file = local_absolute_path(config_file, pwd);
config_dir = fileparts(config_file);
config = jsondecode(fileread(config_file));
local_require_fields(config, {'qe_path', 'epsilon', 'sigma'});

% This function resides at the repository root, so startup no longer
% depends on where the submission script or JSON file resides.
repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);
KSSOLV_startup;

if isfield(config, 'rng_seed') && ~isempty(config.rng_seed)
    rng(config.rng_seed, 'twister');
end

read_vxc = true;
if isfield(config, 'read_vxc')
    read_vxc = logical(config.read_vxc);
end
qe_path = local_absolute_path(char(config.qe_path), config_dir);
if ~isfolder(qe_path)
    error('gw_run_json:QEPathNotFound', 'QE input directory not found: %s', qe_path);
end

fprintf('GW JSON configuration: %s\n', config_file);
fprintf('QE input directory:     %s\n', qe_path);
[sys, options, syms] = read_qe_gw(qe_path, read_vxc);
[sys, options] = gw_setup(sys, options);

eps = config.epsilon;
eps = local_complete_epsilon_bands(eps, options);
eps = epsilon(sys, options, syms, eps);

sig = config.sigma;
sig = sigma(eps, sig, sys, options, syms);

qp_file = fullfile(config_dir, 'qp.dat');
if isfield(config, 'qp_file') && ~isempty(config.qp_file)
    qp_file = local_absolute_path(char(config.qp_file), config_dir);
end
local_write_qp_dat(sig, qp_file);
fprintf('Quasiparticle levels:   %s\n', qp_file);

if isfield(config, 'save_file') && ~isempty(config.save_file)
    save_file = local_absolute_path(char(config.save_file), config_dir);
    save_dir = fileparts(save_file);
    if ~isempty(save_dir) && ~isfolder(save_dir)
        mkdir(save_dir);
    end
    save(save_file, 'eps', 'sig', 'sys', 'options', 'syms', 'config', '-v7.3');
    fprintf('Saved GW results:       %s\n', save_file);
end
end

function eps = local_complete_epsilon_bands(eps, options)
if ~isfield(eps, 'nbnd') || isempty(eps.nbnd)
    error('gw_run_json:MissingEpsilonBands', 'JSON field epsilon.nbnd is required.');
end
if ~isfield(eps, 'nv') || isempty(eps.nv)
    eps.nv = options.nv;
end
if ~isfield(eps, 'nc') || isempty(eps.nc)
    eps.nc = eps.nbnd - eps.nv;
end
if eps.nc <= 0
    error('gw_run_json:InvalidEpsilonBands', ...
        'epsilon.nbnd (%d) must exceed epsilon.nv (%d).', eps.nbnd, eps.nv);
end
end

function local_require_fields(config, fields)
for ii = 1:numel(fields)
    if ~isfield(config, fields{ii}) || isempty(config.(fields{ii}))
        error('gw_run_json:MissingField', 'JSON field "%s" is required.', fields{ii});
    end
end
end

function path_out = local_absolute_path(path_in, base_dir)
if ispc && ~isempty(regexp(path_in, '^[A-Za-z]:[\\/]', 'once'))
    path_out = path_in;
elseif startsWith(path_in, filesep)
    path_out = path_in;
else
    path_out = fullfile(base_dir, path_in);
end
end

function local_write_qp_dat(sig, filename)
% Write quasiparticle levels in BerkeleyGW-style qp.dat form (all energies eV).
required = {'emf', 'ax', 'asx', 'ach', 'achx', 'sig', 'vxc', 'eqp0'};
for ii = 1:numel(required)
    if ~isfield(sig, required{ii})
        error('gw_run_json:MissingQPField', ...
            'sig.%s is required; run sigma before writing qp.dat.', required{ii});
    end
end

filename = char(filename);
output_dir = fileparts(filename);
if ~isempty(output_dir) && ~isfolder(output_dir)
    mkdir(output_dir);
end
fid = fopen(filename, 'w');
if fid < 0
    error('gw_run_json:QPFileOpenFailed', 'Cannot open %s for writing.', filename);
end
cleanup = onCleanup(@() fclose(fid));

has_eqp1 = isfield(sig, 'eqp1') && ~isempty(sig.eqp1);
nspin = size(sig.eqp0, 3);
nkn = size(sig.eqp0, 2);
ndiag = size(sig.eqp0, 1);
if isfield(sig, 'ndiag_min')
    bands = sig.ndiag_min + (0:ndiag - 1);
else
    bands = 1:ndiag;
end
include_k_spin = nkn > 1 || nspin > 1;
if include_k_spin
    fprintf(fid, '%4s %4s %5s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n', ...
        'n', 'k', 'spin', 'Emf', 'Eo', 'X', 'SX-X', 'CH', 'Sig', 'Vxc', 'Eqp0', 'Eqp1');
else
    fprintf(fid, '%4s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n', ...
        'n', 'Emf', 'Eo', 'X', 'SX-X', 'CH', 'Sig', 'Vxc', 'Eqp0', 'Eqp1');
end

for ispin = 1:nspin
    for ik = 1:nkn
        for in = 1:ndiag
            eqp1 = sig.eqp0(in, ik, ispin);
            if has_eqp1
                eqp1 = sig.eqp1(in, ik, ispin);
            end
            emf = sig.emf(in, ik, ispin);
            ch = sig.ach(in, ik, ispin) + sig.achx(in, ik, ispin);
            values = real([emf, emf, sig.ax(in, ik, ispin), ...
                sig.asx(in, ik, ispin), ch, sig.sig(in, ik, ispin), ...
                sig.vxc(in, ik, ispin), sig.eqp0(in, ik, ispin), eqp1]);
            if include_k_spin
                fprintf(fid, '%4d %4d %5d', bands(in), ik, ispin);
            else
                fprintf(fid, '%4d', bands(in));
            end
            fprintf(fid, ' %12.6f', values);
            fprintf(fid, '\n');
        end
    end
end
end
