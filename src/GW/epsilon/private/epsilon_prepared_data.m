function prepared = epsilon_prepared_data( ...
    ctx, iq, ik, wfnk_all, wfnkq_all, fft_all, idx_all)
%EPSILON_PREPARED_DATA Return precomputed epsilon block data when enabled.

if ~ctx.precompute_wav
    prepared = [];
    return;
end

prepared.wfnk = wfnk_all{iq, ik};
prepared.wfnkq = wfnkq_all{iq, ik};
prepared.fft = fft_all{iq};
prepared.idx.k = idx_all.k{iq, ik};
prepared.idx.q = idx_all.q{iq};
prepared.idx.kq = idx_all.kq{iq, ik};
end
