function value = screened_contract(kernel, coeff_left, coeff_right)
%ISDF.SCREENED_CONTRACT Contract a reduced screened Coulomb kernel.

if nargin < 3 || isempty(coeff_right)
    coeff_right = coeff_left;
end
value = coeff_left(:).' * kernel * conj(coeff_right(:));
end
