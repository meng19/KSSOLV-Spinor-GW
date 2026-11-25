function [eps, sig, sys, options, syms] = gw_run_json(config_file)
%GW_RUN_JSON Run a GW calculation described by a portable JSON file.
%   GW_RUN_JSON(CONFIG_FILE) may be called from any working directory.  All
%   relative paths in CONFIG_FILE are interpreted relative to that file.
%
%   Required JSON fields are:
%     qe_path, epsilon, sigma
%   Optional fields are:
%     read_vxc (true), rng_seed, qp_file, save_file,
%     epsilon.checkpoint_file, sigma.epsilon_checkpoint_file
%
%   Strings in CONFIG_FILE may reference environment variables with
%   ${NAME} or %NAME%; referencing an unset variable is an error.  A
%   leading "~" in a path field is expanded to the home directory.
%   Together with the relative-path rule above this lets one JSON file run
%   unchanged on any machine, for example
%     "qe_path": "${QE_DATA_ROOT}/si8"
%
%   Example launchers (the launcher, the JSON, and the repository may live
%   in unrelated directories):
%     run_gw.sh  case.json
%     run_gw.ps1 case.json
%     matlab -batch "addpath(getenv('KSSOLV_ROOT')); \
%         gw_run_json(getenv('GW_JSON'))"

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
config = local_expand_environment(config);
local_require_fields(config, {'qe_path', 'epsilon', 'sigma'});

% Checkpoint files are consumed by epsilon/sigma rather than here, so they
% are resolved now to keep their meaning independent of the working
% directory that happens to run this function.
config.epsilon = local_resolve_optional_paths(config.epsilon, config_dir, ...
    {'checkpoint_file', 'checkpoint', 'breakpoint'}, true);
config.sigma = local_resolve_optional_paths(config.sigma, config_dir, ...
    {'epsilon_checkpoint_file', 'epsilon_checkpoint'}, false);

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

% The output phases are timed separately from epsilon/sigma, so the timers
% accumulated by those stages are cleared before the writes below.
gw_timer('reset');

qp_file = fullfile(config_dir, 'qp.dat');
if isfield(config, 'qp_file') && ~isempty(config.qp_file)
    qp_file = local_absolute_path(char(config.qp_file), config_dir);
end
gw_timer('start', 'Quasiparticle table write');
local_write_qp_dat(sig, qp_file);
gw_timer('stop', 'Quasiparticle table write');
fprintf('Quasiparticle levels:   %s\n', qp_file);

if isfield(config, 'save_file') && ~isempty(config.save_file)
    save_file = local_absolute_path(char(config.save_file), config_dir);
    save_dir = fileparts(save_file);
    if ~isempty(save_dir) && ~isfolder(save_dir)
        mkdir(save_dir);
    end
    gw_timer('start', 'GW results save');
    save(save_file, 'eps', 'sig', 'sys', 'options', 'syms', 'config', '-v7.3');
    gw_timer('stop', 'GW results save');
    fprintf('Saved GW results:       %s\n', save_file);
end
gw_timer('report', 'GW output writing');
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

function container = local_resolve_optional_paths(container, base_dir, fields, create_dir)
% Resolve optional file paths stored inside the epsilon/sigma objects so
% they follow the same rule as qe_path, qp_file, and save_file.
if ~isstruct(container)
    return;
end
for ii = 1:numel(fields)
    name = fields{ii};
    if ~isfield(container, name)
        continue;
    end
    value = container.(name);
    % sig.epsilon_checkpoint may hold a restored eps struct; leave it alone.
    is_text = ischar(value) || (isstring(value) && isscalar(value));
    if ~is_text || isempty(value)
        continue;
    end
    value = local_absolute_path(char(value), base_dir);
    value_dir = fileparts(value);
    if create_dir && ~isempty(value_dir) && ~isfolder(value_dir)
        mkdir(value_dir);
    end
    container.(name) = value;
end
end

function value = local_expand_environment(value)
% Recursively replace ${NAME} and %NAME% with environment variables.
switch class(value)
    case 'struct'
        names = fieldnames(value);
        for ii = 1:numel(names)
            value.(names{ii}) = local_expand_environment(value.(names{ii}));
        end
    case 'cell'
        for ii = 1:numel(value)
            value{ii} = local_expand_environment(value{ii});
        end
    case {'char', 'string'}
        value = local_expand_string(value);
end
end

function value = local_expand_string(value)
% Replace environment references one at a time so that a variable whose
% value itself contains a reference cannot loop forever.
pattern = '\$\{([A-Za-z_][A-Za-z0-9_]*)\}|%([A-Za-z_][A-Za-z0-9_]*)%';
text = char(value);
for iter = 1:1000
    [match, start_idx] = regexp(text, pattern, 'match', 'start', 'once');
    if isempty(match)
        break;
    end
    if match(1) == '$'
        name = match(3:end - 1);
    else
        name = match(2:end - 1);
    end
    replacement = getenv(name);
    if isempty(replacement)
        error('gw_run_json:UndefinedEnvironmentVariable', ...
            ['Environment variable "%s" referenced by the JSON ' ...
             'configuration is not set.'], name);
    end
    text = [text(1:start_idx - 1), replacement, ...
        text(start_idx + numel(match):end)];
end
value = text;
end

function path_out = local_absolute_path(path_in, base_dir)
% Expand a leading "~" before the relative/absolute decision below.
if strcmp(path_in, '~') || startsWith(path_in, ['~' filesep]) || ...
        startsWith(path_in, '~/')
    home_dir = getenv('HOME');
    if isempty(home_dir)
        home_dir = getenv('USERPROFILE');
    end
    if isempty(home_dir)
        error('gw_run_json:MissingHome', ...
            'Cannot expand "~" because HOME and USERPROFILE are unset.');
    end
    path_in = [home_dir path_in(2:end)];
end
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
