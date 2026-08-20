function contribution = sigma_make_contribution( ...
    ctx, asx_loc, ax_loc, ach_loc, achx_loc, ...
    omega, iw_lda, achx_loc_nn)
%SIGMA_MAKE_CONTRIBUTION Pack common sigma contraction outputs.

contribution.asx = asx_loc;
contribution.ax = ax_loc;
contribution.ach = ach_loc;
contribution.achx = achx_loc;
contribution.omega = omega;
contribution.iw_lda = iw_lda;
if ctx.sig.freq_dep == 2
    contribution.asx_freq = asx_loc;
    contribution.ach_freq = ach_loc;
else
    contribution.asx_freq = [];
    contribution.ach_freq = [];
end
contribution.achx_nn = achx_loc_nn;
end
