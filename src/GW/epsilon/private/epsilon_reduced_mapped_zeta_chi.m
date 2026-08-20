function zeta_chi = epsilon_reduced_mapped_zeta_chi(zeta_g, g_maps)
%EPSILON_REDUCED_MAPPED_ZETA_CHI Build mapped zeta blocks for screened W.

zeta_cells = cell(1, numel(g_maps));
for it = 1:numel(g_maps)
    zeta_cells{it} = conj(zeta_g(g_maps{it}, :));
end
zeta_chi = cat(2, zeta_cells{:});
end
