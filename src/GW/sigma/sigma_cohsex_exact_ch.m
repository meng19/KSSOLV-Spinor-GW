function achx_loc = sigma_cohsex_exact_ch(in, ispin, fbz, indrk, iq, aqsch, eps_inv_I_coul, sig, igpp_tmp, valid_indices_tmp)
if sig.freq_dep == 0 || sig.exact_static_ch == 1
    ncouls = fbz.nmtx_cutoff(1, indrk(iq));
    aqsch_tmp = zeros(ncouls);
    aqsch_tmp(valid_indices_tmp) = aqsch{in, ispin}(igpp_tmp(valid_indices_tmp));
    aqsch_tmp(~valid_indices_tmp) = 0;
    achx_loc = aqsch_tmp .* eps_inv_I_coul(:,:,1);
end