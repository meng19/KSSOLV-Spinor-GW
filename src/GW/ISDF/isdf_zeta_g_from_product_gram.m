function [zeta_g, info] = isdf_zeta_g_from_product_gram(c1, c2, idx_q, fftgrid, options)
%ISDF_ZETA_G_FROM_PRODUCT_GRAM Solve the ISDF basis from Gram factors.

[zeta_real, info] = isdf_stable_right_solve(c1, c2, options);
zeta_g = isdf_zeta_real_to_g(zeta_real, idx_q, fftgrid);
end
