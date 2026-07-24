function coulg = coulG_select(coul_options, nmtx, isrtx, ekin, job_type, g_q, gvec, sys, iq)
%COULG_SELECT Select and build the truncated Coulomb potential.

coul_cut = coul_options.coul_cut;
if strcmp(coul_cut, 'spherical_truncation')
    if ~isfield(coul_options, 'coul_cutoff')
        error('coul_cutoff is required when coul_cut is spherical_truncation.');
    end
    coulg = coulG_spherical_truncation(nmtx, isrtx, ekin, ...
        coul_options.coul_cutoff, job_type);
elseif strcmp(coul_cut, 'cell_box_truncation')
    if iq > 1
        error('cell_box_truncation only support one Gamma=0 calculation')
    end
    coulg = coulG_cell_box_truncation(g_q, gvec, sys);
else
    error('Unknown truncation schemes for the Coulomb potential: %s. Please choose spherical_truncation or cell_box_truncation.', coul_cut);
end
end
