function sig = quasi_energy(nspin, ndiag_min, ndiag_max, emf, vxc, sig)
%% Symmetrize matrix elements
sig = symmetrize_sig(sig, emf, ndiag_min, ndiag_max, nspin, {'cor', 'sig'});

for ik = 1:sig.nkn
    %% Get quasiparticle energy
    [~, ind] = ismember(sig.qpt(ik, :), vxc.kpoints, 'rows');
    if (nspin == 2)
        sig.vxc(:, ik, 1) = vxc.value(2*ndiag_min - 1:2:2*ndiag_max, ind);
        sig.vxc(:, ik, 2) = vxc.value(2*ndiag_min:2:2*ndiag_max, ind);
    else
        sig.vxc(:, ik) = vxc.value(ndiag_min:ndiag_max, ind);
    end
end
sig.eqp0 = emf(ndiag_min:ndiag_max, :, :) - sig.vxc + sig.sig;
end