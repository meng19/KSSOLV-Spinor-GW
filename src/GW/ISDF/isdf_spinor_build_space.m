function space = isdf_spinor_build_space(left_components, right_components, idx_q, fftgrid, options)
%ISDF_SPINOR_BUILD_SPACE Build ISDF space for summed spinor products.

products = isdf_spinor_products_from_real(left_components, right_components);
space = isdf_product_space_from_products(products, idx_q, fftgrid, options);
space.left_mu_components = cell(size(left_components));
space.right_mu_components = cell(size(right_components));
for ispinor = 1:numel(left_components)
    space.left_mu_components{ispinor} = left_components{ispinor}(space.ind_mu, :);
    space.right_mu_components{ispinor} = right_components{ispinor}(space.ind_mu, :);
end
end
