function contribution = sigma_make_contribution( ...
    ctx, asx_loc, ax_loc, ach_loc, achx_loc, varargin)
%SIGMA_MAKE_CONTRIBUTION Pack common sigma contraction outputs.

contribution.asx = asx_loc;
contribution.ax = ax_loc;
contribution.ach = ach_loc;
contribution.achx = achx_loc;
if ctx.sig.freq_dep == 2
    omega = varargin{1};
    iw_lda = varargin{2};
    achx_loc_nn = varargin{3};
    contribution.omega = omega;
    contribution.iw_lda = iw_lda;
    contribution.asx_freq = asx_loc;
    contribution.ach_freq = ach_loc;
    contribution.achx_nn = achx_loc_nn;
end
end
