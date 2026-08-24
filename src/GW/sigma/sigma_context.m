function ctx = sigma_context(eps, sig, sys, options, syms, metadata_only)
%SIGMA_CONTEXT Build shared setup data for all sigma methods.

if nargin < 6
    metadata_only = false;
end

ctx.eps = eps;
ctx.sig = sig;
ctx.sys = sys;
ctx.options = options;
ctx.syms = syms;
ctx.method = gw_resolve_method(sig.isdf, 'sigma');
ctx.ryd = 13.6056923;
ctx.nbands = sig.nbnd;
ctx.band_range = sig.ndiag_min:sig.ndiag_max;
ctx.nspin = sys.nspin;
ctx.nspinor = sys.nspinor;
ctx.nk = sys.nkpts;
ctx.wfc_cutoff = 2 * sys.ecut;
ctx.use_gpu = sig.use_gpu && exist('gpuDevice', 'file');
ctx.precompute_wav = sig.precompute_wav;

ctx.sig.qpt = options.kpts;
ctx.sig.nkn = sys.nkpts;
if sig.freq_dep == 2 && sig.freq_dep_method == 2
    if sig.freq_grid_shift == 2
        ctx.sig.nfreq_grid = 2 * fix( ...
            sig.max_freq_eval / sig.delta_freq_eval) + 1;
        ctx.sig.freq_grid = 0:sig.delta_freq_eval:2 * sig.max_freq_eval;
    end
    ctx.sig.nfreq_integral = eps.nfreq;
    ctx.sig.freq_integral = eps.freq;
    ctx.sig.nfreq_integral_imag = eps.nfreq_imag;
    ctx.sig.nfreq_integral_real = eps.nfreq - eps.nfreq_imag;
elseif sig.freq_dep == 0
    ctx.sig.nfreq_grid = 1;
end

sigrid = Ggrid(sys, 4 * sys.ecut);
ctx.gvec = Gvector(sigrid, sys);
ctx.gr = fullbz(options, syms, true);
ctx.fact = 1 / (ctx.gr.nf * sys.vol);
ctx.coulfact = 8 * pi * ctx.fact;
ctx.grid_size = [sys.n1, sys.n2, sys.n3];

if ctx.use_gpu
    try
        ctx.fft.Nfft1 = gpuArray.zeros(ctx.grid_size);
        ctx.fft.Nfft2 = gpuArray.zeros(ctx.grid_size);
    catch
        warning('GPU memory insufficient for FFT grids. Falling back to CPU.');
        ctx.use_gpu = false;
        ctx.fft.Nfft1 = zeros(ctx.grid_size);
        ctx.fft.Nfft2 = zeros(ctx.grid_size);
    end
else
    ctx.fft.Nfft1 = zeros(ctx.grid_size);
    ctx.fft.Nfft2 = zeros(ctx.grid_size);
end
ctx.fft.size = prod(ctx.grid_size);

ctx.ekin_ir = zeros(ctx.gvec.ng, ctx.nk);
for ik = 1:ctx.nk
    rk = ctx.sig.qpt(ik, :);
    [ctx.ekin_ir(:, ik), ctx.sig.isrtx(:, ik)] = sortrx( ...
        rk, ctx.gvec.ng, ctx.gvec.mill, sys);
    ctx.sig.nmtx(:, ik) = gcutoff(ctx.gvec.ng, ...
        ctx.ekin_ir(:, ik), ctx.sig.isrtx(:, ik), eps.cutoff);
    ctx.sig.mtx{:, ik} = ctx.gvec.mill( ...
        ctx.sig.isrtx(1:ctx.sig.nmtx(ik), ik), :);
end

ctx.fbz = struct();
ctx.ekin_fbz = zeros(ctx.gvec.ng, ctx.gr.nf);
ctx.eps_inv_fbz = cell(ctx.gr.nf, 1);
ctx.screened_fbz = cell(ctx.gr.nf, 1);
has_full_inverse = isfield(eps, 'inv') && ~isempty(eps.inv);
has_reduced_w = isfield(eps, 'isdf_screened_w') && ...
    ~isempty(eps.isdf_screened_w);
for iq = 1:ctx.gr.nf
    qq = ctx.gr.f(iq, :);
    [ctx.ekin_fbz(:, iq), ctx.fbz.isrtx(:, iq)] = sortrx( ...
        qq, ctx.gvec.ng, ctx.gvec.mill, sys);
    ctx.fbz.nmtx(:, iq) = gcutoff(ctx.gvec.ng, ...
        ctx.ekin_fbz(:, iq), ctx.fbz.isrtx(:, iq), ctx.wfc_cutoff);
    ctx.fbz.mtx{:, iq} = ctx.gvec.mill( ...
        ctx.fbz.isrtx(1:ctx.fbz.nmtx(iq), iq), :);
    ctx.fbz.nmtx_cutoff(:, iq) = gcutoff(ctx.gvec.ng, ...
        ctx.ekin_fbz(:, iq), ctx.fbz.isrtx(:, iq), eps.cutoff);
    ctx.fbz.mtx_cutoff{:, iq} = ctx.gvec.mill( ...
        ctx.fbz.isrtx(1:ctx.fbz.nmtx_cutoff(iq), iq), :);

    isorti = zeros(ctx.gvec.ng, 1);
    for ig = 1:ctx.gvec.ng
        isorti(ctx.fbz.isrtx(ig, iq)) = ig;
    end
    ctx.fbz.isorti(:, iq) = isorti;

    irq = ctx.gr.indr(iq);
    itran = ctx.gr.itran(iq);
    qk = ctx.gr.r(irq, :) * syms.mtrx{itran, :};
    [~, kgq] = krange(qk, 1e-9);
    isorti_ir = zeros(ctx.gvec.ng, 1);
    for ig = 1:ctx.gvec.ng
        isorti_ir(ctx.sig.isrtx(ig, irq)) = ig;
    end
    indt = gmap(ctx.gvec, syms, ctx.sig.nmtx(irq), itran, kgq, ...
        ctx.fbz.isrtx(:, iq), isorti_ir, sys);
    ctx.fbz.screen_map{iq} = indt;

    if ~metadata_only
        if ~strcmp(ctx.method, 'reduced_basis')
            ctx.eps_inv_fbz{iq} = eps.inv{irq}(indt, indt, :);
        elseif isfield(eps, 'inv') && numel(eps.inv) >= irq && ...
                ~isempty(eps.inv{irq})
            ctx.eps_inv_fbz{iq} = eps.inv{irq}(indt, indt, :);
        end
        if isfield(eps, 'isdf_screened_w') && ...
                numel(eps.isdf_screened_w) >= irq && ...
                ~isempty(eps.isdf_screened_w{irq})
            ctx.screened_fbz{iq} = sigma_map_screened_w( ...
                eps.isdf_screened_w{irq}, indt);
        end
    end
end

if ~metadata_only && strcmp(ctx.method, 'reduced_basis') && ...
        ~has_full_inverse && ~has_reduced_w
    error('ISDF:ReducedSigmaMissingScreening', ...
        'Provide eps.inv or eps.isdf_screened_w for reduced-basis sigma.');
end

ctx.kdata = cell(ctx.nk, 1);
for ik = 1:ctx.nk
    ctx.kdata{ik} = sigma_kdata(ctx, ik);
end
end
