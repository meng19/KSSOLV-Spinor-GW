function eps = isdf_epsilon_reduced_basis(sys, options, syms, eps)
%ISDF_EPSILON_REDUCED_BASIS Static epsilon in a reduced ISDF basis.

if eps.freq_dep ~= 0
    error('ISDF:ReducedEpsilonFrequency', ...
        'ISDF reduced-basis epsilon requires eps.freq_dep = 0.');
end
if eps.use_gpu
    error('ISDF:ReducedEpsilonGPU', ...
        'ISDF reduced-basis epsilon currently supports CPU execution only.');
end

nbands = eps.nbnd;
nspin = sys.nspin;
nspinor = sys.nspinor;
wfc_cutoff = 2 * sys.ecut;

fprintf(['System parameters: nvbands = %d, ncbands = %d, nbands = %d, ' ...
    'nspin = %d, nspinor = %d\n'], ...
    eps.nv, eps.nc, nbands, nspin, nspinor);
fprintf('Initializing static ISDF reduced-basis epsilon calculation\n');

output_mode = lower(eps.isdf.output);
need_full_inverse = any(strcmp(output_mode, {'full_inverse', 'both'}));
need_screened_w = any(strcmp(output_mode, {'screened_w', 'both'}));
if ~need_full_inverse && ~need_screened_w
    error('ISDF:ReducedEpsilonOutput', ...
        'Unknown ISDF reduced-basis epsilon output "%s".', eps.isdf.output);
end

sigrid = Ggrid(sys, 4 * sys.ecut);
gvec = Gvector(sigrid, sys);
gr = fullbz(options, syms, true);

pol.qpt = options.kpts;
pol.nfreq = 1;
pol.freq = 0;
ekin = zeros(gvec.ng, sys.nkpts);
for iq = 1:sys.nkpts
    qq = pol.qpt(iq, :);
    [ekin(:, iq), pol.isrtx(:, iq)] = sortrx( ...
        qq, gvec.ng, gvec.mill, sys);
    pol.nmtx(:, iq) = gcutoff( ...
        gvec.ng, ekin(:, iq), pol.isrtx(:, iq), eps.cutoff);
    pol.mtx{:, iq} = gvec.mill( ...
        pol.isrtx(1:pol.nmtx(iq), iq), :);

    box_min = zeros(1, 3);
    box_max = zeros(1, 3);
    [box_min, box_max] = get_gvecs_bounds( ...
        pol.mtx{:, iq}, box_min, box_max);
    pol.fftgrid{:, iq} = min( ...
        options.wfn_fftgrid + box_max - box_min, options.fftgrid);
end

fact = 4 / (gr.nf * sys.vol * nspin * nspinor);
if need_full_inverse
    eps_inv = cell(sys.nkpts, 1);
end
if need_screened_w
    eps.isdf_screened_w = cell(sys.nkpts, 1);
end

solver_options.method = eps.isdf.reduced_solver;
solver_options.froErr = eps.isdf.cauchy_froErr;
solver_options.MaxIter = eps.isdf.cauchy_MaxIter;

for iq = 1:sys.nkpts
    qq = pol.qpt(iq, :);
    syms_qq = subgrp(qq, syms);
    [nrq, neq, indrk] = irrbz(syms_qq, gr);
    g_maps = local_star_g_maps( ...
        iq, nrq, neq, indrk, syms_qq, syms, gr, gvec, pol, sys);

    nmtx = pol.nmtx(iq);
    if need_full_inverse
        chi0_sum = zeros(nmtx, nmtx);
    end
    zeta_blocks = {};
    coeff_blocks = {};

    fprintf(['\n[Epsilon reduced] Q-point %2d/%2d | ' ...
        'Q-vector = (%8.4f, %8.4f, %8.4f) | nmtx = %d'], ...
        iq, sys.nkpts, qq(1), qq(2), qq(3), nmtx);

    for ispin = 1:nspin
        for ik = 1:nrq
            rk = gr.f(indrk(ik), :);
            wfnk = genwf(rk, gr, gvec, syms, sys, options, ...
                wfc_cutoff, nbands, false);
            wfnkq = genwf(rk + qq, gr, gvec, syms, sys, options, ...
                wfc_cutoff, nbands, false);
            [fft, idx] = epsilon_prefft( ...
                wfnkq, wfnk, iq, ik, pol, [], [], false);

            occ_vkq = get_occ(options, wfnkq.ikq, ispin);
            no_v = sum(occ_vkq > 0);
            occ_ck = get_occ(options, wfnk.ikq, ispin);
            conduction_bands = (sum(occ_ck > 0) + 1):nbands;
            if no_v == 0 || isempty(conduction_bands)
                continue;
            end

            left = cell(1, nspinor);
            right = cell(1, nspinor);
            for ispinor = 1:nspinor
                left{ispinor} = isdf_wavefunction_real_component( ...
                    wfnkq, fft.Nfft1, idx.kq, ispin, ispinor, 1:no_v);
                right{ispinor} = isdf_wavefunction_real_component( ...
                    wfnk, fft.Nfft2, idx.k, ispin, ispinor, ...
                    conduction_bands);
            end

            isdf_options = eps.isdf;
            if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
                npairs = no_v * numel(conduction_bands);
                isdf_options.rank = ceil(sqrt(npairs) * eps.isdf.rank_ratio);
            end
            space = isdf_build_space( ...
                left, right, idx.q, size(fft.Nfft1), isdf_options);

            ev_occ = options.ev(1:no_v, wfnkq.ikq, ispin);
            ev_unocc = options.ev(conduction_bands, wfnk.ikq, ispin);
            polar = isdf_reduced_polarizability( ...
                space, ev_occ, ev_unocc, solver_options);
            eps.isdf_reduced_info{iq, ispin, ik} = polar.info;
            eps.isdf_reduced_rank{iq, ispin, ik} = space.rank;

            % Symmetry-equivalent k points have the same reduced
            % coefficient, but their G rows must be explicitly mapped.
            for it = 1:numel(g_maps{ik})
                zeta_star = space.zeta_g(g_maps{ik}{it}, :);
                zeta_chi = conj(zeta_star);
                coeff_chi = conj(polar.coeff);
                if need_full_inverse
                    chi0_sum = chi0_sum + ...
                        zeta_chi * coeff_chi * zeta_chi';
                end
                if need_screened_w
                    zeta_blocks{end + 1} = zeta_chi; %#ok<AGROW>
                    coeff_blocks{end + 1} = coeff_chi; %#ok<AGROW>
                end
            end
        end
    end

    coulg = coulG_select(eps, nmtx, pol.isrtx(:, iq), ...
        ekin(:, iq), 0, pol.mtx{:, iq}, gvec, sys, iq);
    if need_full_inverse
        epsilon_matrix = eye(nmtx) - coulg(:) .* (fact * chi0_sum);
        eps_inv{iq} = inv(epsilon_matrix);
    end
    if need_screened_w && ~isempty(zeta_blocks)
        combined_space.zeta_g = cat(2, zeta_blocks{:});
        combined_polar.coeff = fact * blkdiag(coeff_blocks{:});
        eps.isdf_screened_w{iq} = isdf_static_screened_interaction( ...
            combined_space, coulg(:), combined_polar);
    end
end

if need_full_inverse
    eps.inv = eps_inv;
end
eps.mtx = pol.mtx;
eps.nmtx = pol.nmtx;
eps.nfreq = pol.nfreq;
eps.freq = pol.freq;

fprintf('\nCalculation of ISDF reduced-basis epsilon completed successfully.\n');
end

function g_maps = local_star_g_maps(iq, nrq, neq, indrk, syms_q, syms, gr, gvec, pol, sys)
% Build the same star-to-G mappings used by the direct epsilon path.

isorti = zeros(gvec.ng, 1);
for ig = 1:gvec.ng
    isorti(pol.isrtx(ig, iq)) = ig;
end

g_maps = cell(nrq, 1);
for ik = 1:nrq
    rk = gr.f(indrk(ik), :);
    [nstar, indst] = rqstar(syms_q, rk);
    if nstar ~= neq(ik)
        error('ISDF:ReducedEpsilonStar', ...
            'K-point star size does not match its irreducible weight.');
    end
    g_maps{ik} = cell(nstar, 1);
    g_maps{ik}{1} = (1:pol.nmtx(iq)).';
    for it = 2:nstar
        itran = syms_q.indsub(indst(it));
        kgq = -syms_q.kgzero(indst(it), :);
        g_maps{ik}{it} = gmap(gvec, syms, pol.nmtx(iq), ...
            itran, kgq, pol.isrtx(:, iq), isorti, sys);
    end
end
end
