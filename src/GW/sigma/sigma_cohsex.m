function [asx_loc, ach_loc] = sigma_cohsex(asx_loc, ach_loc, occ_kq, aqsn, aqsm, eps_inv_I_coul)
aqs_product = aqsn.* aqsm';
aqs_eps_coul = sum(aqs_product .* eps_inv_I_coul, 'all');
%% Calculate SX - X
if occ_kq > 0
    asx_loc = asx_loc - occ_kq * aqs_eps_coul;
end
%% Calculate CH without exact ch correlation
ach_loc = ach_loc + aqs_eps_coul;
end