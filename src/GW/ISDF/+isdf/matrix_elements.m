function gme = matrix_elements(left, right, idx_q, fftgrid, options)
%ISDF.MATRIX_ELEMENTS Fourier matrix elements for component products.

if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
space = isdf.build_space(left, right, idx_q, fftgrid, options);
nleft = size(left{1}, 2);
nright = size(right{1}, 2);
gme = reshape(space.zeta_g * space.product_mu, ...
    numel(idx_q), nleft, nright);
end
