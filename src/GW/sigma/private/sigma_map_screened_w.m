function mapped = sigma_map_screened_w(screened, indt)
%SIGMA_MAP_SCREENED_W Permute screened W into a full-BZ G ordering.

if size(screened.zeta_g, 1) < max(indt) || ...
        numel(screened.epsilon_vcoul) < max(indt)
    error('ISDF:ReducedSigmaMapSize', ...
        'Screened interaction is too small for the requested G mapping.');
end
mapped = screened;
mapped.zeta_g = screened.zeta_g(indt, :);
mapped.epsilon_vcoul = screened.epsilon_vcoul(indt);
end
