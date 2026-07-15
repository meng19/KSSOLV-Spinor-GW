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
