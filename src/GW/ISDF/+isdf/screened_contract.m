function value = screened_contract(kernel, coeff_left, coeff_right)
if nargin < 3
    coeff_right = [];
end
value = isdf_screened_coulomb_contract( ...
    kernel, coeff_left, coeff_right);
end
