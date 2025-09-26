function vxc = cal_vxc(sys, options)
ryd = 13.6056923; % Ry --> eV
[vhart, vxc_r, exc, rho] = getVhxc(sys, options.rho0);
if (sys.nspin == 2)
    vxc_r{1} = 2 * vxc_r{1}; % Ha --> Ry
    vxc_r{2} = 2 * vxc_r{2}; % Ha --> Ry
else
    vxc_r = 2 * vxc_r;
end
Nfft(1) = sys.n1;
Nfft(2) = sys.n2;
Nfft(3) = sys.n3;
vxc.value = zeros(sys.nspin * sys.nbnd, sys.nkpts);

for ik = 1 : sys.nkpts
    g_k = options.mill{1, ik};
    nmtx = size(g_k, 1);
    bidx_k = g2fft_grid(g_k, Nfft(1), Nfft(2), Nfft(3));
    idx_k = sub2ind([Nfft(1), Nfft(2), Nfft(3)], bidx_k(:,1), bidx_k(:,2), bidx_k(:,3));
    Nfft1 = zeros(Nfft(1), Nfft(2), Nfft(3));
    for ib = 1 : sys.nbnd
        for ispin = 1 : sys.nspin
            for ispinor = 1 : sys.nspinor
                evc = options.X0.wavefuncell{ik, ispin}.psi( 1 + (ispinor - 1) * nmtx : ispinor * nmtx, ib);
                Nfft1(:) = 0;
                Nfft1(idx_k) = evc;
                if (sys.nspin == 2)
                    psic = ifftn(Nfft1) * (Nfft(1) * Nfft(2) * Nfft(3)) .* vxc_r{ispin};
                else
                    psic = ifftn(Nfft1) * (Nfft(1) * Nfft(2) * Nfft(3)) .* vxc_r;
                end
                psic = fftn(psic) / (Nfft(1) * Nfft(2) * Nfft(3));
                hpsi = psic(idx_k);
                vxc.value(ispin + (ib - 1) * sys.nspin, ik) = vxc.value(ispin + (ib - 1) * sys.nspin, ik) + sum(conj(evc).* hpsi) * ryd;
            end
        end
    end
end
vxc.value = real(vxc.value);
vxc.kpoints = round(sys.kpts / sys.bvec^2, 6);
end