function polar = polarizability(space, ev_occ, ev_unocc, options)
%ISDF.POLARIZABILITY Build reduced polarizability coefficients.

if nargin < 4 || isempty(options)
    options = struct();
end
if isfield(space, 'left_mu_components') && ...
        isfield(space, 'right_mu_components')
    left_mu = space.left_mu_components;
    right_mu = space.right_mu_components;
elseif isfield(space, 'phi_mu') && isfield(space, 'psi_mu')
    left_mu = {conj(space.phi_mu)};
    right_mu = {space.psi_mu};
else
    error('ISDF:ReducedPolarizabilitySpace', ...
        'space must contain selected scalar or component wavefunctions.');
end

[coeff, info] = local_comega_cstar( ...
    left_mu, right_mu, ev_occ, ev_unocc, options);
polar = struct('coeff', coeff, 'info', info);
end

% ---- Solver dispatch ----

function [result, info] = local_comega_cstar( ...
    left, right, ev_occ, ev_unocc, options)
if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
if ~isfield(options, 'method') || isempty(options.method)
    options.method = 'cauchy';
end
if ~isfield(options, 'freq') || isempty(options.freq)
    options.freq = 0;
end
ev_occ = ev_occ(:);
ev_unocc = ev_unocc(:);
freq = options.freq(:).';

switch lower(options.method)
    case 'direct'
        result = local_direct(left, right, ev_occ, ev_unocc, freq);
        info = struct('method', 'direct', 'iterations', 0, ...
            'relative_error', 0);
    case 'cauchy'
        if numel(freq) ~= 1 || freq ~= 0
            error('ISDF:CauchyFullFrequencyUnsupported', ...
                ['Cauchy reduced polarizability currently supports ' ...
                 'static frequency only. Use reduced_solver = direct.']);
        end
        if ~isfield(options, 'froErr')
            options.froErr = 1e-8;
        end
        if ~isfield(options, 'MaxIter')
            options.MaxIter = 12;
        end
        [result, relative_error, iterations] = local_cauchy( ...
            left, right, ev_occ, ev_unocc, options);
        info = struct('method', 'cauchy', ...
            'iterations', iterations, 'relative_error', relative_error);
    otherwise
        error('ISDF:UnknownCOmegaMethod', ...
            'Unknown COmegaCstar method "%s".', options.method);
end
end

% ---- Direct frequency-page polarizability ----

function result = local_direct(left, right, ev_occ, ev_unocc, freq)
nmu = size(left{1}, 1);
nv = numel(ev_occ);
nc = numel(ev_unocc);
ncomponents = numel(left);
if ncomponents ~= numel(right)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end
result = zeros(nmu, nmu, numel(freq));
for iv = 1:nv
    for ic = 1:nc
        coefficient = zeros(nmu, 1);
        for icomponent = 1:ncomponents
            coefficient = coefficient + ...
                conj(left{icomponent}(:, iv)) .* right{icomponent}(:, ic);
        end
        energy_diff = ev_occ(iv) - ev_unocc(ic);
        for ifreq = 1:numel(freq)
            if freq(ifreq) == 0
                eden = 1 / energy_diff;
            else
                eden = energy_diff / (energy_diff^2 - freq(ifreq)^2);
            end
            result(:, :, ifreq) = result(:, :, ifreq) + ...
                coefficient * coefficient' * eden;
        end
    end
end
end

% ---- Static Cauchy polarizability ----

function [result, relative_error, iteration] = local_cauchy( ...
    left, right, ev_occ, ev_unocc, options)
gap_min = ev_unocc(1) - ev_occ(end);
if gap_min <= 0
    error('ISDF:CauchyNoGap', ...
        'Cauchy integral requires a positive band gap.');
end
center = 0.5 * (ev_occ(1) + ev_occ(end));
half_width = 0.5 * (ev_occ(end) - ev_occ(1));
radius = half_width + 0.5 * gap_min;
max_radius = ev_unocc(1) - center;
if radius >= max_radius
    radius = 0.9 * max_radius;
end

nmu = size(left{1}, 1);
ncomponents = numel(left);
if ncomponents ~= numel(right)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end
previous = [];
relative_error = inf;
for iteration = 1:options.MaxIter
    npoints = 2^(iteration + 3);
    result = zeros(nmu, nmu);
    for ipoint = 0:npoints-1
        theta = 2 * pi * ipoint / npoints;
        exp_theta = exp(1i * theta);
        z = center + radius * exp_theta;
        for icomponent = 1:ncomponents
            left_s = left{icomponent};
            right_s = right{icomponent};
            for jcomponent = 1:ncomponents
                left_t = left{jcomponent};
                right_t = right{jcomponent};
                occ_matrix = conj(left_s) * ...
                    diag(1 ./ (z - ev_occ)) * left_t.';
                unocc_matrix = right_s * ...
                    diag(1 ./ (z - ev_unocc)) * right_t';
                result = result + (occ_matrix .* unocc_matrix) * ...
                    (radius * exp_theta / npoints);
            end
        end
    end
    if ~isempty(previous)
        relative_error = norm(result - previous, 'fro') / ...
            max(1, norm(result, 'fro'));
        if relative_error <= options.froErr
            return;
        end
    end
    previous = result;
end
end
