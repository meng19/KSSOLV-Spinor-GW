function kernel = screened_kernel(screened, target_zeta_g, contract_vcoul)
kernel = isdf_screened_coulomb_kernel( ...
    screened, target_zeta_g, contract_vcoul);
end
