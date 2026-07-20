function eps = isdf_epsilon_cauchy_polarizability(sys, options, syms, eps)
%ISDF_EPSILON_CAUCHY_POLARIZABILITY Static epsilon via component ISDF/Cauchy.

if eps.freq_dep ~= 0
    error('ISDF:CauchyEpsilonFrequency', ...
        'ISDF Cauchy epsilon path requires eps.freq_dep = 0.');
end

if eps.use_gpu
    error('ISDF:CauchyEpsilonGPU', ...
        'ISDF Cauchy epsilon path currently supports CPU execution only.');
end

if sys.nkpts ~= 1
    error('ISDF:CauchyEpsilonKPoints', ...
        'Initial ISDF Cauchy epsilon implementation supports single-k systems only.');
end

ryd = 13.6056923;
nvbands = eps.nv;
ncbands = eps.nc;
nbands = eps.nbnd;
nspin = sys.nspin;
nspinor = sys.nspinor;
wfc_cutoff = sys.ecut * 2;

fprintf('System parameters: nvbands = %d, ncbands = %d, nbands = %d, nspin = %d, nspinor = %d\n', ...
    nvbands, ncbands, nbands, nspin, nspinor);
fprintf('Initializing static ISDF Cauchy epsilon calculation (frequency = 0)\n');

sigrid = Ggrid(sys, 4 * sys.ecut);
gvec = Gvector(sigrid, sys);
pol.qpt = options.kpts;
pol.nfreq = 1;
pol.freq = 0;

gr = fullbz(options, syms, true);
ekin = zeros(gvec.ng, sys.nkpts);
for iq = 1:sys.nkpts
    qq = pol.qpt(iq, :);
    [ekin(:, iq), pol.isrtx(:, iq)] = sortrx(qq, gvec.ng, gvec.mill, sys);
    pol.nmtx(:, iq) = gcutoff(gvec.ng, ekin(:, iq), pol.isrtx(:, iq), eps.cutoff);
    pol.mtx{:, iq} = gvec.mill(pol.isrtx(1:pol.nmtx(iq), iq), :);

    eps_box_min = zeros([1 3]);
    eps_box_max = zeros([1 3]);
    [eps_box_min, eps_box_max] = get_gvecs_bounds(pol.mtx{:, iq}, eps_box_min, eps_box_max);
    pol.fftgrid{:, iq} = min((options.wfn_fftgrid + eps_box_max - eps_box_min), options.fftgrid);
end

store_full_inverse = eps.isdf.store_full_inverse;
if store_full_inverse
    eps_inv = cell(sys.nkpts, 1);
end
fact = 4 / (gr.nf * sys.vol * nspin * nspinor);

for iq = 1:sys.nkpts
    qq = pol.qpt(iq, :);
    if norm(qq) > 1e-12
        error('ISDF:CauchyEpsilonQPoint', ...
            'Initial ISDF Cauchy epsilon implementation supports Gamma q only.');
    end

    syms_qq = subgrp(qq, syms);
    [nrq, neq, indrk] = irrbz(syms_qq, gr);
    if nrq ~= 1 || neq(1) ~= 1
        error('ISDF:CauchyEpsilonSymmetry', ...
            'Initial ISDF Cauchy epsilon implementation supports one irreducible k point only.');
    end

    nmtx_current = pol.nmtx(iq);
    chi0_sum = zeros(nmtx_current, nmtx_current);
    screened_current = [];

    fprintf('\n[Epsilon Cauchy] K-point %2d/%2d | K-vector = (%8.4f, %8.4f, %8.4f) | nmtx = %2d', ...
        iq, sys.nkpts, qq(1), qq(2), qq(3), nmtx_current);

    for ispin = 1:nspin
        rk = gr.f(indrk(1), :);
        wfnk = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff, nbands, false);
        rkq = rk + qq;
        wfnkq = genwf(rkq, gr, gvec, syms, sys, options, wfc_cutoff, nbands, false);
        [fft, idx] = epsilon_prefft(wfnkq, wfnk, iq, 1, pol, [], [], false);

        occ_vkq = get_occ(options, wfnkq.ikq, ispin);
        no_v = sum(occ_vkq > 0);
        occ_ck = get_occ(options, wfnk.ikq, ispin);
        no_c_start = sum(occ_ck > 0) + 1;
        conduction_bands = no_c_start:nbands;

        if no_v == 0 || isempty(conduction_bands)
            continue;
        end

        left_components = cell(1, nspinor);
        right_components = cell(1, nspinor);
        for ispinor = 1:nspinor
            left_components{ispinor} = isdf_wavefunction_real_component( ...
                wfnkq, fft.Nfft1, idx.kq, ispin, ispinor, 1:no_v);
            right_components{ispinor} = isdf_wavefunction_real_component( ...
                wfnk, fft.Nfft2, idx.k, ispin, ispinor, conduction_bands);
        end

        isdf_options = eps.isdf;
        if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
            isdf_options.rank = ceil(sqrt(no_v * numel(conduction_bands)) * eps.isdf.rank_ratio);
        end
        vc_space = isdf_build_space(left_components, right_components, ...
            idx.q, size(fft.Nfft1), isdf_options);

        cauchy_options.method = eps.isdf.cauchy_method;
        cauchy_options.froErr = eps.isdf.cauchy_froErr;
        cauchy_options.MaxIter = eps.isdf.cauchy_MaxIter;
        ev_occ = options.ev(1:no_v, wfnkq.ikq, ispin);
        ev_unocc = options.ev(conduction_bands, wfnk.ikq, ispin);
        [coeff, info] = isdf_comega_cstar(vc_space.left_mu_components, ...
            vc_space.right_mu_components, ev_occ, ev_unocc, cauchy_options);

        chi0_sum = chi0_sum + conj(vc_space.zeta_g) * conj(coeff) * vc_space.zeta_g.';
        eps.isdf_cauchy_info{iq, ispin} = info;
        eps.isdf_cauchy_rank{iq, ispin} = vc_space.rank;
        if nspin == 1
            screened_space = vc_space;
            screened_space.zeta_g = conj(vc_space.zeta_g);
            screened_polar = struct();
            screened_polar.coeff = conj(coeff) * fact;
            screened_current = struct('space', screened_space, 'polar', screened_polar);
        end
    end

    chi0_sum = chi0_sum * fact;

    if strcmp(eps.coul_cut, 'spherical_truncation')
        coulg = coulG_spherical_truncation(nmtx_current, pol.isrtx(:, iq), ...
            ekin(:, iq), eps.coul_cutoff, 0);
    elseif strcmp(eps.coul_cut, 'cell_box_truncation')
        coulg = coulG_cell_box_truncation(pol.mtx{:, iq}, gvec, sys);
    else
        error('ISDF:CauchyEpsilonCoulomb', ...
            'Unknown truncation scheme "%s".', eps.coul_cut);
    end

    if store_full_inverse
        eps_tmp = eye(nmtx_current) - bsxfun(@times, coulg(:), chi0_sum);
        eps_inv{iq} = inv(eps_tmp);
    end
    if ~isempty(screened_current)
        eps.isdf_screened{iq, 1} = isdf_static_screened_interaction( ...
            screened_current.space, coulg(:), screened_current.polar);
    end
end

if store_full_inverse
    eps.inv = eps_inv;
end
eps.mtx = pol.mtx;
eps.nmtx = pol.nmtx;
eps.nfreq = pol.nfreq;
eps.freq = pol.freq;

fprintf('\nCalculation of ISDF Cauchy epsilon completed successfully.\n');
end
