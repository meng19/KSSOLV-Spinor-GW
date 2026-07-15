function value = isdf_screened_coulomb_contract(kernel, coeff_left, coeff_right)
%ISDF_SCREENED_COULOMB_CONTRACT Contract reduced screened Coulomb kernel.

if nargin < 3 || isempty(coeff_right)
    coeff_right = coeff_left;
end

value = coeff_left(:).' * kernel * conj(coeff_right(:));
end
