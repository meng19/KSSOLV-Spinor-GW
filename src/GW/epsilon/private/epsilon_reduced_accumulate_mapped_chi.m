function chi0 = epsilon_reduced_accumulate_mapped_chi( ...
    chi0, zeta_g, coeff_chi, g_maps)
%EPSILON_REDUCED_ACCUMULATE_MAPPED_CHI Accumulate reduced chi0 pages.

base_zeta = conj(zeta_g);
for ifreq = 1:size(coeff_chi, 3)
    base_chi = base_zeta * coeff_chi(:, :, ifreq) * base_zeta';
    for it = 1:numel(g_maps)
        gmap = g_maps{it};
        chi0(:, :, ifreq) = chi0(:, :, ifreq) + base_chi(gmap, gmap);
    end
end
end
