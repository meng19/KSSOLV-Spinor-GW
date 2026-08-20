function acc = epsilon_reduced_accumulate(~, acc, contribution, block)
%EPSILON_REDUCED_ACCUMULATE Accumulate reduced-basis block outputs.

coeff_chi = conj(contribution.polar.coeff);
if acc.need_full_inverse
    acc.chi0 = epsilon_reduced_accumulate_mapped_chi( ...
        acc.chi0, contribution.space.zeta_g, coeff_chi, block.g_maps);
end
if acc.need_screened_w
    zeta_chi = epsilon_reduced_mapped_zeta_chi( ...
        contribution.space.zeta_g, block.g_maps);
    acc.zeta_blocks{end + 1} = zeta_chi;
    acc.coeff_blocks{end + 1} = epsilon_repeat_page_blkdiag( ...
        coeff_chi, numel(block.g_maps));
end
acc.rank{block.ispin, block.ik} = contribution.space.rank;
acc.info{block.ispin, block.ik} = contribution.polar.info;
end
