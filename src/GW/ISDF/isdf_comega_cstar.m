function [result, info] = isdf_comega_cstar(phi, psi, ev_occ, ev_unocc, options)
%ISDF_COMEGA_CSTAR Compute C*Omega^{-1}*C' in ISDF interpolation space.
%   C_mu,vc = conj(phi_mu,v) * psi_mu,c and
%   Omega_vc = ev_occ(v) - ev_unocc(c).

if nargin < 5 || isempty(options)
    options = struct();
end

if ~iscell(phi)
    phi = {phi};
end
if ~iscell(psi)
    psi = {psi};
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

function result = isdf_comega_cstar_direct(left_mu_components, right_mu_components, ev_occ, ev_unocc)
nmu = size(left_mu_components{1}, 1);
nv = numel(ev_occ);
nc = numel(ev_unocc);
ncomponents = numel(left_mu_components);

if ncomponents ~= numel(right_mu_components)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end

result = zeros(nmu, nmu);
for iv = 1:nv
    for ic = 1:nc
        c = zeros(nmu, 1);
        for icomponent = 1:ncomponents
            c = c + conj(left_mu_components{icomponent}(:, iv)) .* ...
                right_mu_components{icomponent}(:, ic);
        end
        result = result + (c * c') / (ev_occ(iv) - ev_unocc(ic));
    end
end
end

function [result, rel_error, iter] = isdf_comega_cstar_cauchy(left_mu_components, right_mu_components, ev_occ, ev_unocc, options)
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

nmu = size(left_mu_components{1}, 1);
ncomponents = numel(left_mu_components);
if ncomponents ~= numel(right_mu_components)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end

previous = [];
rel_error = inf;

for iter = 1:options.MaxIter
    npts = 2^(iter + 3);
    result = zeros(nmu, nmu);
    for ipt = 0:npts-1
        theta = 2 * pi * ipt / npts;
        exp_theta = exp(1i * theta);
        z = center + radius * exp_theta;
        for icomponent = 1:ncomponents
            left_s = left_mu_components{icomponent};
            right_s = right_mu_components{icomponent};
            for jcomponent = 1:ncomponents
                left_t = left_mu_components{jcomponent};
                right_t = right_mu_components{jcomponent};
                occ_matrix = conj(left_s) * diag(1 ./ (z - ev_occ)) * left_t.';
                unocc_matrix = right_s * diag(1 ./ (z - ev_unocc)) * right_t';
                result = result + (occ_matrix .* unocc_matrix) * ...
                    (radius * exp_theta / npts);
            end
        end
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
