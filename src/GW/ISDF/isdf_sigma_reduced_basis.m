function sig = isdf_sigma_reduced_basis(eps, sig, sys, options, syms)
%ISDF_SIGMA_REDUCED_BASIS Static COHSEX self-energy in an ISDF basis.

if eps.freq_dep ~= 0 || sig.freq_dep ~= 0
    error('ISDF:ReducedSigmaFrequency', ...
        'ISDF reduced-basis sigma requires static epsilon and sigma.');
end
if isfield(sig, 'use_gpu') && sig.use_gpu
    error('ISDF:ReducedSigmaGPU', ...
        'ISDF reduced-basis sigma currently supports CPU execution only.');
end

ryd = 13.6056923;
nbands = sig.nbnd;
ndiag_min = sig.ndiag_min;
ndiag_max = sig.ndiag_max;
ndiag = ndiag_max - ndiag_min + 1;
nspin = sys.nspin;
nspinor = sys.nspinor;
wfc_cutoff = 2 * sys.ecut;

sigrid = Ggrid(sys, 4 * sys.ecut);
gvec = Gvector(sigrid, sys);
gr = fullbz(options, syms, true);
fact = 1 / (gr.nf * sys.vol);

sig.qpt = options.kpts;
sig.nkn = sys.nkpts;
ekin_ir = zeros(gvec.ng, sys.nkpts);
for ik = 1:sig.nkn
    rk = sig.qpt(ik, :);
    [ekin_ir(:, ik), sig.isrtx(:, ik)] = sortrx( ...
        rk, gvec.ng, gvec.mill, sys);
    sig.nmtx(:, ik) = gcutoff( ...
        gvec.ng, ekin_ir(:, ik), sig.isrtx(:, ik), eps.cutoff);
    sig.mtx{:, ik} = gvec.mill( ...
        sig.isrtx(1:sig.nmtx(ik), ik), :);
end

ekin_fbz = zeros(gvec.ng, gr.nf);
for iq = 1:gr.nf
    qq = gr.f(iq, :);
    [ekin_fbz(:, iq), fbz.isrtx(:, iq)] = sortrx( ...
        qq, gvec.ng, gvec.mill, sys);
    fbz.nmtx(:, iq) = gcutoff( ...
        gvec.ng, ekin_fbz(:, iq), fbz.isrtx(:, iq), wfc_cutoff);
    fbz.mtx{:, iq} = gvec.mill( ...
        fbz.isrtx(1:fbz.nmtx(iq), iq), :);
    fbz.nmtx_cutoff(:, iq) = gcutoff( ...
        gvec.ng, ekin_fbz(:, iq), fbz.isrtx(:, iq), eps.cutoff);
    fbz.mtx_cutoff{:, iq} = gvec.mill( ...
        fbz.isrtx(1:fbz.nmtx_cutoff(iq), iq), :);
    isorti = zeros(gvec.ng, 1);
    for ig = 1:gvec.ng
        isorti(fbz.isrtx(ig, iq)) = ig;
    end
    fbz.isorti(:, iq) = isorti;
end

has_full_inverse = isfield(eps, 'inv') && ~isempty(eps.inv);
has_reduced_w = isfield(eps, 'isdf_screened_w') && ...
    ~isempty(eps.isdf_screened_w);
if ~has_full_inverse && ~has_reduced_w
    error('ISDF:ReducedSigmaMissingScreening', ...
        'Provide eps.inv or eps.isdf_screened_w for reduced-basis sigma.');
end

eps_inv_fbz = cell(gr.nf, 1);
screened_fbz = cell(gr.nf, 1);
for iq = 1:gr.nf
    irq = gr.indr(iq);
    itran = gr.itran(iq);
    qk = gr.r(irq, :) * syms.mtrx{itran, :};
    [~, kgq] = krange(qk, 1e-9);
    isorti = zeros(gvec.ng, 1);
    for ig = 1:gvec.ng
        isorti(sig.isrtx(ig, irq)) = ig;
    end
    indt = gmap(gvec, syms, sig.nmtx(irq), itran, kgq, ...
        fbz.isrtx(:, iq), isorti, sys);

    if has_full_inverse && numel(eps.inv) >= irq && ~isempty(eps.inv{irq})
        eps_inv_fbz{iq} = eps.inv{irq}(indt, indt, :);
    end
    if has_reduced_w && numel(eps.isdf_screened_w) >= irq && ...
            ~isempty(eps.isdf_screened_w{irq})
        screened_fbz{iq} = local_map_screened_w( ...
            eps.isdf_screened_w{irq}, indt);
    end
end

grid_size = [sys.n1, sys.n2, sys.n3];
fft.Nfft1 = zeros(grid_size);
fft.Nfft2 = zeros(grid_size);
fft.size = prod(grid_size);

asx = zeros([ndiag, sys.nkpts, nspin]);
ax = zeros([ndiag, sys.nkpts, nspin]);
ach = zeros([ndiag, sys.nkpts, nspin]);
achx = zeros([ndiag, sys.nkpts, nspin]);
use_exact_static_ch = isfield(sig, 'exact_static_ch') && sig.exact_static_ch;
no_q_symmetry = isfield(sig, 'no_symmetries_q_grid') && ...
    sig.no_symmetries_q_grid;

fprintf('Starting ISDF reduced-basis static sigma calculation...\n');
for ispin = 1:nspin
    for in = ndiag_min:ndiag_max
        n_index = in - ndiag_min + 1;
        for ik = 1:sig.nkn
            rk = sig.qpt(ik, :);
            syms_rk = subgrp(rk, syms);
            [nrk, neq, indrk] = irrbz(syms_rk, gr);
            if no_q_symmetry
                nrk = gr.nf;
                indrk = 1:nrk;
                neq = ones(1, nrk);
            end

            wfnk = genwf(rk, gr, gvec, syms, sys, options, ...
                wfc_cutoff, nbands, false);
            asx_sum = 0;
            ax_sum = 0;
            ach_sum = 0;
            achx_sum = 0;
            aqsch = cell(nbands, nspin);

            for iq = 1:nrk
                iq_fbz = indrk(iq);
                qq = gr.f(iq_fbz, :);
                if ~no_q_symmetry
                    [nstar, ~] = rqstar(syms_rk, qq);
                    if nstar ~= neq(iq)
                        error('ISDF:ReducedSigmaStar', ...
                            'Q-point star size does not match its weight.');
                    end
                end

                n_cutoff = fbz.nmtx_cutoff(iq_fbz);
                coulg = coulG_select(sig, fbz.nmtx(iq_fbz), ...
                    fbz.isrtx(:, iq_fbz), ekin_fbz(:, iq_fbz), 1, ...
                    fbz.mtx{:, iq_fbz}, gvec, sys, iq_fbz);
                coulg_cutoff = coulg(1:n_cutoff);

                wfnkq = genwf(rk - qq, gr, gvec, syms, sys, options, ...
                    wfc_cutoff, nbands, false);
                idx = sigma_prefft(wfnkq, wfnk, fbz.mtx{:, iq_fbz}, ...
                    iq, ik, sys, [], false);
                occ_kq = get_occ(options, wfnkq.ikq, ispin);

                left = cell(1, nspinor);
                right = cell(1, nspinor);
                for ispinor = 1:nspinor
                    left{ispinor} = isdf_wavefunction_real_component( ...
                        wfnk, fft.Nfft1, idx.k, ispin, ispinor, in);
                    right{ispinor} = isdf_wavefunction_real_component( ...
                        wfnkq, fft.Nfft2, idx.kq, ispin, ispinor, 1:nbands);
                end

                isdf_options = sig.isdf;
                if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
                    isdf_options.rank = ceil( ...
                        sqrt(nbands) * sig.isdf.rank_ratio);
                end
                space = isdf_build_space( ...
                    left, right, idx.q, grid_size, isdf_options);
                aqs_isdf = space.zeta_g * space.product_mu;

                use_reduced_w = ~isempty(screened_fbz{iq_fbz});
                if use_reduced_w
                    screened = screened_fbz{iq_fbz};
                    if size(screened.zeta_g, 1) ~= n_cutoff
                        error('ISDF:ReducedSigmaScreenedSize', ...
                            'Reduced screened interaction does not match sigma cutoff.');
                    end
                    target_zeta = space.zeta_g(1:n_cutoff, :);
                    screened_kernel = isdf_screened_coulomb_kernel( ...
                        screened, target_zeta, coulg_cutoff);
                    if use_exact_static_ch
                        screened_matrix = fact * ...
                            isdf_screened_coulomb_kernel( ...
                            screened, [], coulg_cutoff);
                    end
                else
                    if isempty(eps_inv_fbz{iq_fbz})
                        error('ISDF:ReducedSigmaMissingQPoint', ...
                            'No screened interaction is available for full-BZ q-point %d.', ...
                            iq_fbz);
                    end
                    eps_inv = eps_inv_fbz{iq_fbz};
                    screened_matrix = fact * ((eps_inv - eye(n_cutoff)) .* ...
                        coulg_cutoff.');
                end

                asx_loc = 0;
                ax_loc = 0;
                ach_loc = 0;
                for nn = 1:nbands
                    aqs = aqs_isdf(:, nn);
                    if occ_kq(nn) > 0
                        ax_loc = ax_loc - occ_kq(nn) * fact * ...
                            sum(abs(aqs).^2 .* coulg);
                    end
                    if use_reduced_w
                        coeff = space.product_mu(:, nn);
                        screened_value = fact * ...
                            isdf_screened_coulomb_contract( ...
                            screened_kernel, coeff);
                        if occ_kq(nn) > 0
                            asx_loc = asx_loc - occ_kq(nn) * screened_value;
                        end
                        ach_loc = ach_loc + screened_value;
                    else
                        aqs_cutoff = aqs(1:n_cutoff);
                        [asx_loc, ach_loc] = sigma_cohsex( ...
                            asx_loc, ach_loc, occ_kq(nn), ...
                            aqs_cutoff, aqs_cutoff, screened_matrix);
                    end
                end

                asx_sum = asx_sum + neq(iq) * asx_loc;
                ax_sum = ax_sum + neq(iq) * ax_loc;
                ach_sum = ach_sum + neq(iq) * ach_loc;

                if use_exact_static_ch
                    if iq_fbz == 1
                        aqsch{in, ispin} = aqs_isdf(:, in);
                    end
                    [igpp, valid_indices] = pre_exact_static_ch( ...
                        fbz, gvec, indrk, iq, false);
                    achx_loc = sigma_cohsex_exact_ch( ...
                        in, ispin, fbz, indrk, iq, aqsch, ...
                        screened_matrix, sig, igpp, valid_indices);
                    achx_sum = achx_sum + neq(iq) * sum(achx_loc, 'all');
                end
            end

            asx(n_index, ik, ispin) = asx_sum;
            ax(n_index, ik, ispin) = ax_sum;
            ach(n_index, ik, ispin) = 0.5 * ach_sum;
            if use_exact_static_ch
                achx(n_index, ik, ispin) = 0.5 * achx_sum;
            end
        end
    end
end

if use_exact_static_ch
    sig.cor = real(asx + achx) * ryd;
    sig.sig = real(asx + ax + achx) * ryd;
else
    sig.cor = real(asx + ach) * ryd;
    sig.sig = real(asx + ax + ach) * ryd;
end
emf = ryd * options.ev;
sig = quasi_energy(nspin, ndiag_min, ndiag_max, emf, sys.vxc, sig);
fprintf('ISDF reduced-basis static sigma calculation completed.\n');
end

function mapped = local_map_screened_w(screened, indt)
% Permute an irreducible-q screened operator into a full-BZ G ordering.

if size(screened.zeta_g, 1) < max(indt) || ...
        numel(screened.epsilon_vcoul) < max(indt)
    error('ISDF:ReducedSigmaMapSize', ...
        'Screened interaction is too small for the requested G mapping.');
end
mapped = screened;
mapped.zeta_g = screened.zeta_g(indt, :);
mapped.epsilon_vcoul = screened.epsilon_vcoul(indt);
end
