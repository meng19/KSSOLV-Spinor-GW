function sig = get_eqp1(nspin, ndiag_min, ndiag_max, emf, omega, iw_lda, asx_freq, ach_freq, achx, ax, sig)
ryd = 13.6056923;
for ispin = 1:nspin
    for ik = 1:sig.nkn
        for in = ndiag_min:ndiag_max
            n_index = in - ndiag_min + 1;
            freqs = squeeze(omega(in, ik, ispin, :));
            iw = iw_lda(in, ik, ispin);
            ax_tmp = ax(n_index,ik,ispin) * ryd;
            asx_tmp = asx_freq{n_index,ik,ispin} * ryd;
            ach_tmp = ach_freq{n_index,ik,ispin} * ryd;
            achx_tmp = achx(n_index,ik,ispin) * ryd;
            asigt_cor = ach_tmp + achx_tmp;
            fmin = emf(in) + ax_tmp + asx_tmp + ach_tmp + achx_tmp - sig.vxc(n_index, ik, ispin) - freqs;
            
            [eqp1, neqp1] = cal_eqp1(freqs, fmin, sig.eqp0(n_index,ik,ispin));
            sig.eqp1(n_index, ik, ispin) = eqp1;
            sig.neqp1(n_index, ik, ispin) = neqp1;
        end
    end
end
end

function [eqp1, neqp1] = cal_eqp1(freqs, fmin, eqp0)
% GET_EQP1 Optimized version - fast vectorized root finding for quasiparticle energy
%   Finds roots of real(fmin(w)) by detecting + to - crossings with linear interp.
%   If multiple roots, keeps the one closest to eqp0 (real part).
%   If none, tries linear extrapolation left/right, else falls back to eqp0.

    nf = length(freqs);
    if nf < 2
        eqp1  = eqp0;
        neqp1 = 0;
        return;
    end

    rf = real(fmin);
    imf = imag(fmin);
    dw = freqs(2) - freqs(1);   % assume uniform grid (as in original)

    % Vectorized detection of sign changes: + -> -
    sign_change = (rf(1:end-1) > 0) & (rf(2:end) < 0);
    idx = find(sign_change);                    % indices i where crossing between i and i+1

    nsols = length(idx);

    if nsols > 0
        % Linear interpolation (real + imag) for all crossings at once
        denom = rf(idx) - rf(idx+1);
        rsoln = freqs(idx) + rf(idx) .* dw ./ denom;

        imsol = ((freqs(idx+1) - rsoln) .* imf(idx) + ...
                 (rsoln - freqs(idx))     .* imf(idx+1)) / dw;

        solns = complex(rsoln, imsol);

        % Select the one closest to eqp0 (real part distance)
        [~, isol] = min(abs(real(solns) - real(eqp0)));
        eqp1  = solns(isol);
        neqp1 = nsols;

    else
        % No crossing → try extrapolation (left or right)
        if rf(1) < 0 && rf(1) > rf(2)          % left extrapolation
            rsoln = (rf(2)*freqs(1) - rf(1)*freqs(2)) / (rf(2) - rf(1));
            imsol = ((freqs(2) - rsoln)*imf(1) + (rsoln - freqs(1))*imf(2)) / dw;
            eqp1  = complex(rsoln, imsol);
            neqp1 = -1;

        elseif rf(nf) > 0 && rf(nf-1) > rf(nf) % right extrapolation
            rsoln = (rf(nf-1)*freqs(nf) - rf(nf)*freqs(nf-1)) / (rf(nf-1) - rf(nf));
            imsol = ((freqs(nf) - rsoln)*imf(nf-1) + (rsoln - freqs(nf-1))*imf(nf)) / dw;
            eqp1  = complex(rsoln, imsol);
            neqp1 = -2;

        else
            eqp1  = eqp0;
            neqp1 = 0;
        end
    end
end