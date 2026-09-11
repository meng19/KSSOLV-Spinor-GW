function vxc = getvxc(sys, options, syms)

[vhart, vxc, exc, rho] = getVhxc(sys, rho);
sigrid = Ggrid(sys, sys.ecut2);
gvec = Gvector(sigrid, sys);


wfc_cutoff = sys.ecut * 2; % Ha --> Ry
gr = fullbz(sys, syms, true);
[ekin, isrtc_kq] = sortrx(rkq, gvec.ng, gvec.mill, sys);
nmtx = gcutoff(gvec.ng, ekin, isrtc_kq, wfc_cutoff);
wfnkq.mill = gvec.mill(isrtc_kq(1:nmtx), :);
for ik = 1 : gr.nf
    for ispin = 1 : sys.nspin
    for ispinor = 1 : sys.nspinor
        wfn_ispinor = X0.wavefuncell{ikrkq, ispin}.psi( 1 + (ispinor - 1) * nmtx : ispinor * nmtx, :);
    end
end
    g_k = wfnk.mill;
    bidx_k = g2fft_grid(g_k, Nfft(1), Nfft(2), Nfft(3));
    idx_k = sub2ind([Nfft(1), Nfft(2), Nfft(3)], bidx_k(:,1), bidx_k(:,2), bidx_k(:,3));
    Nfft1(idx_k) = wfnk.psi{ispinor};
    Nfft2 = fftn(conj(fftn(Nfft1)) .* fftn(Nfft2)) / (Nfft(1)*Nfft(2)*Nfft(3));
end
end