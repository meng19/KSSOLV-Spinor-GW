function chi0 = epsilon_add_state_batch(chi0, gme, eden, g_maps)
%EPSILON_ADD_STATE_BATCH Accumulate one matrix-element state batch.
%
% In save_mem mode g_maps is provided, so star members are accumulated as
% G-index permutations without expanding gme. In deferred mode gme is
% already star-expanded, so g_maps is empty and one larger GEMM is added
% per frequency page.

if nargin < 4
    g_maps = [];
end
gme = reshape(gme, size(gme, 1), []);
eden = reshape(eden, [], size(eden, 3));
has_maps = ~isempty(g_maps);
for ifreq = 1:size(eden, 2)
    weight = eden(:, ifreq).';
    if isa(gme, 'gpuArray') && ~isa(weight, 'gpuArray')
        weight = gpuArray(weight);
    end
    weighted_gme = bsxfun(@times, gme, weight);
    base_chi = conj(gme) * weighted_gme.';
    if has_maps
        for it = 1:numel(g_maps)
            gmap = g_maps{it};
            chi0(:, :, ifreq) = chi0(:, :, ifreq) + ...
                base_chi(gmap, gmap);
        end
    else
        chi0(:, :, ifreq) = chi0(:, :, ifreq) + base_chi;
    end
end
end
