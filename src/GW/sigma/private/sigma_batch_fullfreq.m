function [asx_loc, ach_loc, achx_row, omega, iw_lda] = sigma_batch_fullfreq( ...
        ctx, block, coeff_all, kernel, progress_work)
%SIGMA_BATCH_FULLFREQ Batch full-frequency self-energy contraction.
%   Evaluates the screened-exchange residue and the Coulomb-hole
%   imaginary-axis contributions for every summation band at once.  Each
%   band's rank-by-rank outer product is folded into block GEMMs against
%   the frequency-page kernel instead of one band at a time.
%
%   The result matches per-band evaluation of SIGMA_FULLFREQ up to BLAS
%   block-summation order.  Unlike the per-band helper, each residue point
%   uses its own interpolation weights (the per-band helper applied the
%   first point's weights to every point through an implicit-expansion
%   row/column mismatch).

sig = ctx.sig;
ryd = ctx.ryd;
nfreqeval = sig.nfreq_grid;
delta_freq = sig.delta_freq_eval;
nfreq_real = sig.nfreq_integral_real;
nfreq_imag = sig.nfreq_integral_imag;
nfreq = nfreq_real + nfreq_imag;
freq_integral = sig.freq_integral;
emf = ctx.options.ev;
ik = block.wfnk.ikq;
ikq = block.wfnkq.ikq;
ispin = block.ispin;
e_nk = emf(block.in, ik, ispin) * ryd;
e_n1kq = emf(1:ctx.nbands, ikq, ispin) * ryd;
occ = reshape(block.occ_kq, 1, []);

% The real-frequency evaluation grid is independent of the summation band.
freq0 = e_nk - delta_freq * (nfreqeval - 1) / 2;
iw_grid = 0:nfreqeval - 1;
omega = freq0 + iw_grid * delta_freq;
wxi = omega.' - e_n1kq.';          % nfreqeval x nbands
wxi(wxi == 0) = 1e-6;
[~, iw_lda] = min(abs(omega - e_nk));
iw_lda = double(iw_lda);

nbands = ctx.nbands;
rank = size(coeff_all, 1);
asx_loc = zeros(nfreqeval, 1, 'like', kernel);
ach_loc = zeros(nfreqeval, 1, 'like', kernel);
achx_row = zeros(1, nbands, 'like', kernel);
kernel_flat = reshape(kernel, rank * rank, nfreq);
% One batch holds rank^2 x batch_size complex coefficients, bounded by
% sig.isdf.fullfreq_batch_bytes (default 48 GiB).
if isfield(ctx.sig.isdf, 'fullfreq_batch_bytes') && ...
        ~isempty(ctx.sig.isdf.fullfreq_batch_bytes)
    batch_budget = ctx.sig.isdf.fullfreq_batch_bytes;
else
    batch_budget = 48 * (1024^3);
end
batch_size = max(1, min(nbands, floor( ...
    batch_budget / max(1, 16 * rank * rank))));
for first = 1:batch_size:nbands
    last = min(nbands, first + batch_size - 1);
    b = last - first + 1;
    coeff = coeff_all(:, first:last);
    % aqsmn = aqsm .* aqsn' in sigma_fullfreq: element (i,j) is
    % coeff_i * conj(coeff_j).
    products = reshape(coeff, [rank, 1, b]) .* ...
        conj(reshape(coeff, [1, rank, b]));
    products = reshape(products, rank * rank, b);
    % Non-conjugate transpose: sigma_fullfreq contracts M_flat.' with the
    % frequency pages.
    svalues = products.' * kernel_flat;         % b x nfreq
    achx_row(first:last) = -0.25 * svalues(:, 1).';

    % Residue contributions on the real-frequency axis.
    wxi_b = wxi(:, first:last);
    need_res = (wxi_b >= 0) ~= (occ(first:last) > 0);
    [irow, icol] = find(need_res);
    if ~isempty(irow)
        wx_abs = abs(wxi_b(need_res));
        ifreq_page = discretize(wx_abs, freq_integral(1:nfreq_real));
        ifreq_page(isnan(ifreq_page)) = nfreq_real - 1;
        % Column indexing into a row-vector frequency table preserves the
        % row orientation; keep every residue vector as a column.
        f_low = freq_integral(ifreq_page).';
        df = freq_integral(ifreq_page + 1).' - f_low;
        fact2 = (wx_abs - f_low) ./ df;
        fact1 = 1 - fact2;
        sres_need = fact1 .* svalues(sub2ind(size(svalues), icol, ifreq_page)) + ...
            fact2 .* svalues(sub2ind(size(svalues), icol, ifreq_page + 1));
        occ_sign = sign(wxi_b(need_res));
        occ_sign(occ_sign == 0) = 1;
        sres_b = zeros(nfreqeval, b, 'like', kernel);
        sres_b(need_res) = occ_sign .* sres_need;
        asx_loc = asx_loc + sum(sres_b, 2);
    end

    % Imaginary-frequency Coulomb-hole contribution (midpoint trapezoid).
    svalues_imag = -svalues(:, nfreq_real + 1:nfreq_real + nfreq_imag);
    imag_freqs = imag(freq_integral(nfreq_real + 1:nfreq_real + nfreq_imag));
    sint_b = zeros(nfreqeval, b, 'like', kernel);
    fstart = imag_freqs(1);
    fmid = (imag_freqs(1) + imag_freqs(2)) / 2;
    sint_b = sint_b + svalues_imag(:, 1).' .* ...
        (atan(fmid ./ wxi_b) - atan(fstart ./ wxi_b));
    for iq = 2:nfreq_imag - 1
        fstart = (imag_freqs(iq - 1) + imag_freqs(iq)) / 2;
        fmid = (imag_freqs(iq) + imag_freqs(iq + 1)) / 2;
        sint_b = sint_b + svalues_imag(:, iq).' .* ...
            (atan(fmid ./ wxi_b) - atan(fstart ./ wxi_b));
    end
    fstart = (imag_freqs(end - 1) + imag_freqs(end)) / 2;
    fend = (-imag_freqs(end - 1) + 3 * imag_freqs(end)) / 2;
    sint_b = sint_b + svalues_imag(:, end).' .* ...
        (atan(fend ./ wxi_b) - atan(fstart ./ wxi_b));
    ach_loc = ach_loc + sum(sint_b, 2) / pi;

    gw_block_progress(block, progress_work * (0.5 + 0.5 * last / nbands), ...
        sprintf('S b%d i%d q%d n%d/%d', ...
        block.in, block.ik, block.iq, last, nbands));
end
end
