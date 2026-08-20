function mapped = epsilon_mapped_gme(gme, g_maps)
%EPSILON_MAPPED_GME Expand gme over equivalent G-space reorderings.

mapped_cells = cell(1, numel(g_maps));
for it = 1:numel(g_maps)
    mapped_cells{it} = gme(g_maps{it}, :, :);
end
mapped = cat(2, mapped_cells{:});
end
