function gw_write_qp_dat(sig, filename)
%GW_WRITE_QP_DAT Write quasiparticle levels in BerkeleyGW-style qp.dat form.
%   The energy columns are in eV.  Eqp1 equals Eqp0 for static calculations,
%   where no frequency-dependent quasiparticle root is available.

if nargin ~= 2 || ~isstruct(sig) || ~(ischar(filename) || isstring(filename))
    error('gw_write_qp_dat:InvalidInput', 'Usage: gw_write_qp_dat(sig, filename).');
end
required = {'emf', 'ax', 'asx', 'ach', 'achx', 'sig', 'vxc', 'eqp0'};
for ii = 1:numel(required)
    if ~isfield(sig, required{ii})
        error('gw_write_qp_dat:MissingField', ...
            'sig.%s is required; run sigma before writing qp.dat.', required{ii});
    end
end

filename = char(filename);
output_dir = fileparts(filename);
if ~isempty(output_dir) && ~isfolder(output_dir)
    mkdir(output_dir);
end

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

fid = fopen(filename, 'w');
if fid < 0
    error('gw_write_qp_dat:OpenFailed', 'Cannot open %s for writing.', filename);
end
cleanup = onCleanup(@() fclose(fid));
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
