function prepared = sigma_prepared_data( ...
    ctx, ik, iq, wfnk, wfnkq_all, idx_all, igpp, valid_indices)
%SIGMA_PREPARED_DATA Return sigma block data available before preparation.

prepared = struct();
prepared.wfnk = wfnk;
if ~ctx.precompute_wav
    return;
end

prepared.wfnkq = wfnkq_all{iq, ik};
prepared.idx.k = idx_all.k{ik};
prepared.idx.q = idx_all.q{iq, ik};
prepared.idx.kq = idx_all.kq{iq, ik};
prepared.igpp = igpp{iq, ik};
prepared.valid_indices = valid_indices{iq, ik};
end
