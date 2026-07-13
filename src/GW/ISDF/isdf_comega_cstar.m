function [result, info] = isdf_comega_cstar(phi, psi, ev_occ, ev_unocc, options)
%ISDF_COMEGA_CSTAR Compute C*Omega^{-1}*C' in ISDF interpolation space.
%   C_mu,vc = conj(phi_mu,v) * psi_mu,c and
%   Omega_vc = ev_occ(v) - ev_unocc(c).

if nargin < 5 || isempty(options)
    options = struct();
end
if ~isfield(options, 'method') || isempty(options.method)
    options.method = 'cauchy';
end

ev_occ = ev_occ(:);
ev_unocc = ev_unocc(:);

switch lower(options.method)
    case 'direct'
        result = isdf_comega_cstar_direct(phi, psi, ev_occ, ev_unocc);
        info = struct('method', 'direct', 'iterations', 0, 'relative_error', 0);
    case 'cauchy'
        if ~isfield(options, 'froErr')
            options.froErr = 1e-8;
        end
        if ~isfield(options, 'MaxIter')
            options.MaxIter = 12;
        end
        [result, rel_error, iterations] = isdf_comega_cstar_cauchy( ...
            phi, psi, ev_occ, ev_unocc, options);
        info = struct('method', 'cauchy', ...
            'iterations', iterations, 'relative_error', rel_error);
    otherwise
        error('ISDF:UnknownCOmegaMethod', ...
            'Unknown COmegaCstar method "%s".', options.method);
end
end

function result = isdf_comega_cstar_direct(phi, psi, ev_occ, ev_unocc)
nmu = size(phi, 1);
result = zeros(nmu, nmu);
for iv = 1:length(ev_occ)
    for ic = 1:length(ev_unocc)
        c = conj(phi(:, iv)) .* psi(:, ic);
        result = result + (c * c') / (ev_occ(iv) - ev_unocc(ic));
    end
end
end

function [result, rel_error, iter] = isdf_comega_cstar_cauchy(phi, psi, ev_occ, ev_unocc, options)
gap_min = ev_unocc(1) - ev_occ(end);
if gap_min <= 0
    error('ISDF:CauchyNoGap', 'Cauchy integral requires a positive band gap.');
end

center = 0.5 * (ev_occ(1) + ev_occ(end));
half_width = 0.5 * (ev_occ(end) - ev_occ(1));
radius = half_width + 0.5 * gap_min;
max_radius = ev_unocc(1) - center;
if radius >= max_radius
    radius = 0.9 * max_radius;
end

previous = [];
rel_error = inf;

for iter = 1:options.MaxIter
    npts = 2^(iter + 3);
    result = zeros(size(phi, 1), size(phi, 1));
    for ipt = 0:npts-1
        theta = 2 * pi * ipt / npts;
        exp_theta = exp(1i * theta);
        z = center + radius * exp_theta;
        occ_matrix = conj(phi) * diag(1 ./ (z - ev_occ)) * phi.';
        unocc_matrix = psi * diag(1 ./ (z - ev_unocc)) * psi';
        result = result + (occ_matrix .* unocc_matrix) * (radius * exp_theta / npts);
    end

    if ~isempty(previous)
        rel_error = norm(result - previous, 'fro') / max(1, norm(result, 'fro'));
        if rel_error <= options.froErr
            return;
        end
    end
    previous = result;
end
end

function [lambda, dlambda] = isdf_cauchy_integrand(t, k, gap_min, gap_max)
L = -log(k) / pi;
[sn, cn, dn] = isdf_ellipjc(t, L);
lambda = sqrt(gap_min * gap_max) .* ((1 / k + sn) ./ (1 / k - sn));
dlambda = cn .* dn .* sqrt(gap_min * gap_max) .* ...
    ((2 / k) ./ (1 / k - sn).^2);
end

function [K, Kprime] = isdf_ellipk_modulus(k)
if k == 0
    K = pi / 2;
    Kprime = inf;
    return;
end
if k == 1
    K = inf;
    Kprime = pi / 2;
    return;
end

[K, ~] = ellipke(k^2);
[Kprime, ~] = ellipke(1 - k^2);
end

function [sn, cn, dn] = isdf_ellipjc(u, L, flag)
if nargin < 3
    [~, Kp] = isdf_ellipkkp(L);
    high = imag(u) > Kp / 2;
    u(high) = 1i * Kp - u(high);
    m = exp(-2 * pi * L);
else
    high = zeros(size(u));
    m = L;
end

if m < 4 * eps
    sinu = sin(u);
    cosu = cos(u);
    sn = sinu + m / 4 * (sinu .* cosu - u) .* cosu;
    cn = cosu + m / 4 * (-sinu .* cosu + u) .* sinu;
    dn = 1 + m / 4 * (cosu.^2 - sinu.^2 - 1);
else
    if m > 1e-3
        kappa = (1 - sqrt(1 - m)) / (1 + sqrt(1 - m));
    else
        kappa = polyval([132, 42, 14, 5, 2, 1, 0], m / 4);
    end
    mu = kappa^2;
    v = u / (1 + kappa);
    [sn1, cn1, dn1] = isdf_ellipjc(v, mu, 1);
    denom = 1 + kappa * sn1.^2;
    sn = (1 + kappa) * sn1 ./ denom;
    cn = cn1 .* dn1 ./ denom;
    dn = (1 - kappa * sn1.^2) ./ denom;
end

if any(high(:))
    snh = sn(high);
    cnh = cn(high);
    dnh = dn(high);
    sn(high) = -1 ./ (sqrt(m) * snh);
    cn(high) = 1i * dnh ./ (sqrt(m) * snh);
    dn(high) = 1i * cnh ./ snh;
end
end

function [K, Kp] = isdf_ellipkkp(L)
if L > 10
    K = pi / 2;
    Kp = pi * L + log(4);
    return;
end

m = exp(-2 * pi * L);
a0 = 1;
b0 = sqrt(1 - m);
mm = 1;
while mm > eps
    a1 = (a0 + b0) / 2;
    b1 = sqrt(a0 .* b0);
    c1 = (a0 - b0) / 2;
    mm = max(abs(c1(:)));
    a0 = a1;
    b0 = b1;
end
K = pi ./ (2 * a1);

a0 = 1;
b0 = sqrt(m);
mm = 1;
while mm > eps
    a1 = (a0 + b0) / 2;
    b1 = sqrt(a0 .* b0);
    c1 = (a0 - b0) / 2;
    mm = max(abs(c1(:)));
    a0 = a1;
    b0 = b1;
end
Kp = pi ./ (2 * a1);
end
