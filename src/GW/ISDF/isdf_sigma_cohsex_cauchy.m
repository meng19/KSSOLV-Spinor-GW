function sig = isdf_sigma_cohsex_cauchy(eps, sig, sys, options, syms)
%ISDF_SIGMA_COHSEX_CAUCHY Static COHSEX self-energy with ISDF/Cauchy.

if eps.freq_dep ~= 0 || sig.freq_dep ~= 0
    error('ISDF:CauchyCOHSEXFrequency', ...
        'ISDF Cauchy COHSEX requires eps.freq_dep = 0 and sig.freq_dep = 0.');
end

if isfield(sig, 'use_gpu') && sig.use_gpu
    error('ISDF:CauchyCOHSEXGPU', ...
        'ISDF Cauchy COHSEX path currently supports CPU execution only.');
end

if sys.nkpts ~= 1
    error('ISDF:CauchyCOHSEXKPoints', ...
        'Initial ISDF Cauchy COHSEX implementation supports single-k systems only.');
end

if isfield(sig, 'exact_static_ch') && sig.exact_static_ch
    error('ISDF:CauchyCOHSEXExactCH', ...
        'ISDF Cauchy COHSEX path currently supports exact_static_ch = false.');
end

sig = isdf_sigma_cohsex_cauchy_singlek(eps, sig, sys, options, syms);
end

function sig = isdf_sigma_cohsex_cauchy_singlek(eps, sig, sys, options, syms)
ryd = 13.6056923;
nbands = sig.nbnd;
ndiag_min = sig.ndiag_min;
ndiag_max = sig.ndiag_max;
ndiag = ndiag_max - ndiag_min + 1;
nspin = sys.nspin;
nspinor = sys.nspinor;
wfc_cutoff = sys.ecut * 2;

sigrid = Ggrid(sys, 4 * sys.ecut);
gvec = Gvector(sigrid, sys);
gr = fullbz(options, syms, true);

fact = 1 / (gr.nf * sys.vol);
sig.qpt = options.kpts;
sig.nkn = sys.nkpts;

for ik = 1:sig.nkn
    rk = sig.qpt(ik, :);
    [ekin(:, ik), sig.isrtx(:, ik)] = sortrx(rk, gvec.ng, gvec.mill, sys);
    sig.nmtx(:, ik) = gcutoff(gvec.ng, ekin(:, ik), sig.isrtx(:, ik), eps.cutoff);
    sig.mtx{:, ik} = gvec.mill(sig.isrtx(1:sig.nmtx(ik), ik), :);
end

for ik = 1:gr.nf
    rk = gr.f(ik, :);
    [ekin_fbz(:, ik), fbz.isrtx(:, ik)] = sortrx(rk, gvec.ng, gvec.mill, sys);
    fbz.nmtx(:, ik) = gcutoff(gvec.ng, ekin_fbz(:, ik), fbz.isrtx(:, ik), wfc_cutoff);
    fbz.mtx{:, ik} = gvec.mill(fbz.isrtx(1:fbz.nmtx(ik), ik), :);
    fbz.nmtx_cutoff(:, ik) = gcutoff(gvec.ng, ekin_fbz(:, ik), fbz.isrtx(:, ik), eps.cutoff);
    fbz.mtx_cutoff{:, ik} = gvec.mill(fbz.isrtx(1:fbz.nmtx_cutoff(ik), ik), :);
end

grid_size = [sys.n1, sys.n2, sys.n3];
fft.Nfft1 = zeros(grid_size);
fft.Nfft2 = zeros(grid_size);
fft.size = prod(grid_size);

asx = zeros([ndiag sys.nkpts nspin]);
ax = zeros([ndiag sys.nkpts nspin]);
ach = zeros([ndiag sys.nkpts nspin]);

fprintf('Starting ISDF Cauchy COHSEX spinor sigma calculation...\n');

for ispin = 1:nspin
    for in = ndiag_min:ndiag_max
        n_index = in - ndiag_min + 1;
        for ik = 1:sig.nkn
            rk = sig.qpt(ik, :);
            qq = gr.f(1, :);
            wfnk = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff, nbands, false);
            rkq = rk - qq;
            wfnkq = genwf(rkq, gr, gvec, syms, sys, options, wfc_cutoff, nbands, false);
            idx = sigma_prefft(wfnkq, wfnk, fbz.mtx{:, 1}, 1, ik, sys, [], false);

            n_cutoff = fbz.nmtx_cutoff(1, 1);
            if strcmp(sig.coul_cut, 'spherical_truncation')
                coulg = coulG_spherical_truncation(fbz.nmtx(1, 1), ...
                    fbz.isrtx(:, 1), ekin_fbz(:, 1), sig.coul_cutoff, 1);
            elseif strcmp(sig.coul_cut, 'cell_box_truncation')
                coulg = coulG_cell_box_truncation(fbz.mtx{:, 1}, gvec, sys);
            else
                error('ISDF:CauchyCOHSEXCoulomb', ...
                    'Unknown truncation scheme "%s".', sig.coul_cut);
            end

            coulg_nocut = fact * coulg;
            coulg_cutoff = coulg_nocut(1:n_cutoff, 1);
            use_reduced_screened = isfield(eps, 'isdf_screened') && ...
                numel(eps.isdf_screened) >= ispin && ~isempty(eps.isdf_screened{1, ispin});
            if ~use_reduced_screened
                if ~isfield(eps, 'inv') || isempty(eps.inv)
                    error('ISDF:CauchyCOHSEXMissingScreened', ...
                        ['eps.inv is unavailable and eps.isdf_screened is missing. ' ...
                        'Use eps.isdf.store_full_inverse = true or provide reduced screened interaction.']);
                end
                eps_inv = eps.inv{1};
                eps_inv_i_coul = (eps_inv - eye(n_cutoff)) .* coulg_cutoff';
            else
                eps_inv_i_coul = [];
            end

            occ_kq = get_occ(options, wfnkq.ikq, ispin);

            left_components = cell(1, nspinor);
            right_components = cell(1, nspinor);
            for ispinor = 1:nspinor
                left_components{ispinor} = isdf_wavefunction_real_component( ...
                    wfnk, fft.Nfft1, idx.k, ispin, ispinor, in);
                right_components{ispinor} = isdf_wavefunction_real_component( ...
                    wfnkq, fft.Nfft2, idx.kq, ispin, ispinor, 1:nbands);
            end

            isdf_options = sig.isdf;
            if ~isfield(isdf_options, 'rank') || isempty(isdf_options.rank)
                isdf_options.rank = ceil(sqrt(nbands) * sig.isdf.rank_ratio);
            end
            sigma_space = isdf_spinor_build_space(left_components, right_components, ...
                idx.q, grid_size, isdf_options);
            aqs_isdf = sigma_space.zeta_g * sigma_space.product_mu;
            if use_reduced_screened
                screened = eps.isdf_screened{1, ispin};
                target_zeta = sigma_space.zeta_g(1:n_cutoff, :);
                if size(screened.zeta_g, 1) ~= n_cutoff
                    error('ISDF:CauchyCOHSEXScreenedSize', ...
                        'Reduced screened interaction size does not match sigma cutoff.');
                end
                screened_kernel = isdf_screened_coulomb_kernel( ...
                    screened, target_zeta, coulg_cutoff);
            else
                screened_kernel = [];
            end

            [cauchy_coeff, cauchy_info] = local_spinor_vc_cauchy( ...
                wfnkq, wfnk, fft, idx, ispin, nspinor, nbands, ...
                options, sig, grid_size);
            sig.isdf_cauchy_coeff{ispin, ik} = cauchy_coeff;
            sig.isdf_cauchy_info{ispin, ik} = cauchy_info;

            asx_loc = 0;
            ax_loc = 0;
            ach_loc = 0;
            for nn = 1:nbands
                aqs_nocut = aqs_isdf(:, nn);
                if occ_kq(nn) > 0
                    ax_loc = ax_loc - occ_kq(nn) * sum(abs(aqs_nocut).^2 .* coulg_nocut);
                end
                if use_reduced_screened
                    coeff = sigma_space.product_mu(:, nn);
                    aqs_eps_coul = isdf_screened_coulomb_contract(screened_kernel, coeff);
                    if occ_kq(nn) > 0
                        asx_loc = asx_loc - occ_kq(nn) * aqs_eps_coul;
                    end
                    ach_loc = ach_loc + aqs_eps_coul;
                else
                    aqs_cutoff = aqs_nocut(1:n_cutoff, 1);
                    [asx_loc, ach_loc] = sigma_cohsex(asx_loc, ach_loc, ...
                        occ_kq(nn), aqs_cutoff, aqs_cutoff, eps_inv_i_coul);
                end
            end

            asx(n_index, ik, ispin) = asx_loc;
            ax(n_index, ik, ispin) = ax_loc;
            ach(n_index, ik, ispin) = 0.5 * ach_loc;
        end
    end
end

sig.cor = real(asx + ach) * ryd;
sig.sig = real(asx + ax + ach) * ryd;

emf = ryd * options.ev;
sig = quasi_energy(nspin, ndiag_min, ndiag_max, emf, sys.vxc, sig);

fprintf('ISDF Cauchy COHSEX spinor sigma calculation completed.\n');
end

function [coeff, info] = local_spinor_vc_cauchy(wfnkq, wfnk, fft, idx, ispin, nspinor, nbands, options, sig, grid_size)
occ_vkq = get_occ(options, wfnkq.ikq, ispin);
no_v = sum(occ_vkq > 0);
occ_ck = get_occ(options, wfnk.ikq, ispin);
no_c_start = sum(occ_ck > 0) + 1;
conduction_bands = no_c_start:nbands;

if no_v == 0 || isempty(conduction_bands)
    coeff = [];
    info = struct('method', sig.isdf.cauchy_method, ...
        'iterations', 0, 'relative_error', 0);
    return;
end

left_components = cell(1, nspinor);
right_components = cell(1, nspinor);
for ispinor = 1:nspinor
    left_components{ispinor} = isdf_wavefunction_real_component( ...
        wfnkq, fft.Nfft1, idx.kq, ispin, ispinor, 1:no_v);
    right_components{ispinor} = isdf_wavefunction_real_component( ...
        wfnk, fft.Nfft2, idx.k, ispin, ispinor, conduction_bands);
end

vc_options = sig.isdf;
if ~isfield(vc_options, 'rank') || isempty(vc_options.rank)
    vc_options.rank = ceil(sqrt(no_v * numel(conduction_bands)) * sig.isdf.rank_ratio);
end
vc_space = isdf_spinor_build_space(left_components, right_components, ...
    idx.q, grid_size, vc_options);

cauchy_options.method = sig.isdf.cauchy_method;
cauchy_options.froErr = sig.isdf.cauchy_froErr;
cauchy_options.MaxIter = sig.isdf.cauchy_MaxIter;
ev_occ = options.ev(1:no_v, wfnkq.ikq, ispin);
ev_unocc = options.ev(conduction_bands, wfnk.ikq, ispin);
[coeff, info] = isdf_spinor_comega_cstar(vc_space.left_mu_components, ...
    vc_space.right_mu_components, ev_occ, ev_unocc, cauchy_options);
end
