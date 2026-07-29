function [ind_mu, products] = sample_points(left, right, options)
%SAMPLE_POINTS Select ISDF interpolation points for component products.

if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
rng(options.seed, 'twister');
products = [];

switch lower(options.sample_method)
    case 'qrcp'
        products = component_products(left, right, [], []);
        ind_mu = local_qrcp(products, options.rank);
    case {'qrcp_randomized', 'randomized_qrcp', 'default'}
        if numel(left) == 1
            ind_mu = local_scalar_randomized(left{1}, right{1}, options);
        else
            npairs = size(left{1}, 2) * size(right{1}, 2);
            projection_rank = max(options.rank, ...
                ceil(options.random_oversampling * options.rank));
            projection_rank = min(projection_rank, npairs);
            projection = randn(npairs, projection_rank);
            has_complex = any(cellfun(@(x) ~isreal(x), left)) || ...
                any(cellfun(@(x) ~isreal(x), right));
            if has_complex
                projection = projection + 1i * ...
                    randn(npairs, projection_rank);
            end
            compressed_products = component_products( ...
                left, right, [], projection);
            ind_mu = local_qrcp(compressed_products, options.rank);
        end
    case 'kmeans'
        if numel(left) == 1
            weight = local_scalar_weight(left{1}, right{1}, options);
        else
            weight = local_component_weight(left, right);
        end
        ind_mu = local_kmeans(weight, options);
    otherwise
        error('ISDF:UnknownSampleMethod', ...
            ['Unknown ISDF sample_method "%s". Supported methods: qrcp, ' ...
             'qrcp_randomized, kmeans, default.'], options.sample_method);
end
end

% ---- Product weights and randomized compression ----

function weight = local_component_weight(left, right)
% Compute sum(abs(P).^2, 2) without allocating the product matrix P.
ngrid = size(left{1}, 1);
nleft = size(left{1}, 2);
nright = size(right{1}, 2);
weight = zeros(ngrid, 1);
for icomponent = 1:numel(left)
    left_i = left{icomponent};
    right_i = right{icomponent};
    if size(left_i, 1) ~= ngrid || size(right_i, 1) ~= ngrid || ...
            size(left_i, 2) ~= nleft || size(right_i, 2) ~= nright
        error('ISDF:ComponentSizeMismatch', ...
            'All components must share grid and band dimensions.');
    end
    for jcomponent = 1:numel(left)
        left_j = left{jcomponent};
        right_j = right{jcomponent};
        if size(left_j, 1) ~= ngrid || size(right_j, 1) ~= ngrid || ...
                size(left_j, 2) ~= nleft || size(right_j, 2) ~= nright
            error('ISDF:ComponentSizeMismatch', ...
                'All components must share grid and band dimensions.');
        end
        weight = weight + sum(conj(left_i) .* left_j, 2) .* ...
            sum(right_i .* conj(right_j), 2);
    end
end
weight = max(real(weight), 0);
end

function ind_mu = local_scalar_randomized(left, right, options)
nleft = size(left, 2);
nright = size(right, 2);
rank_mu = options.rank;
sample_rank = max(rank_mu, ...
    ceil(options.random_oversampling * rank_mu));
left_rank = min(max(1, ceil(sqrt((nleft / nright) * sample_rank))), ...
    nleft);
right_rank = min(max(1, ceil(sqrt((nright / nleft) * sample_rank))), ...
    nright);

left_projection = randn(nleft, left_rank);
right_projection = randn(nright, right_rank);
if ~isreal(left) || ~isreal(right)
    left_projection = left_projection + 1i * randn(nleft, left_rank);
    right_projection = right_projection + 1i * randn(nright, right_rank);
end
compressed_left = conj(left) * left_projection;
compressed_right = right * right_projection;
products = local_pair_products(compressed_left, compressed_right);
ind_mu = local_qrcp(products, rank_mu);
end

% ---- QRCP sampling ----

function products = local_pair_products(left, right)
[ngrid, nleft] = size(left);
nright = size(right, 2);
products = zeros(ngrid, nleft * nright);
for iright = 1:nright
    for ileft = 1:nleft
        products(:, ileft + (iright - 1) * nleft) = ...
            left(:, ileft) .* right(:, iright);
    end
end
end

function ind_mu = local_qrcp(products, rank_mu)
[~, ~, pivots] = qr(products', 0);
ind_mu = reshape(pivots(1:rank_mu), 1, []);
end

% ---- K-means weights and center selection ----

function weight = local_scalar_weight(left, right, options)
rho_left = sum(abs(left).^2, 2);
rho_right = sum(abs(right).^2, 2);
switch lower(options.weight)
    case {'prod', 'multiply'}
        weight = rho_left .* rho_right;
    case 'add'
        weight = rho_left + rho_right;
    case 'power'
        weight = (rho_left .* rho_right).^(options.power / 2);
    otherwise
        error('ISDF:UnknownWeight', ...
            ['Unknown ISDF weight "%s". Supported weights: ' ...
             'prod, add, power.'], options.weight);
end
weight = real(weight(:));
weight(~isfinite(weight) | weight < 0) = 0;
if all(weight == 0)
    weight = ones(size(weight));
end
end

function ind_mu = local_kmeans(weight, options)
rank_mu = options.rank;
weight = real(weight(:));
points = local_grid_points(numel(weight), options);
ngrid = size(points, 1);
if rank_mu > ngrid
    error('ISDF:RankTooLarge', ...
        'ISDF rank cannot exceed the number of grid points.');
end

switch lower(options.init)
    case 'sequence'
        ind_mu = round(linspace(1, ngrid, rank_mu)).';
    case 'random'
        ind_mu = randperm(ngrid, rank_mu).';
    case 'wrs'
        ind_mu = local_weighted_sample(abs(weight), rank_mu, false);
    case {'kmeans++', 'kmeans++_false'}
        ind_mu = local_kmeanspp_init(points, abs(weight), ...
            rank_mu, options);
    otherwise
        error('ISDF:UnknownKmeansInit', ...
            'Unknown ISDF K-means init "%s".', options.init);
end

ind_mu = local_unique_fill(ind_mu, weight, rank_mu, ngrid);
if options.kmeans_max_iter == 0
    ind_mu = reshape(ind_mu, 1, []);
    return;
end

last_ind = ind_mu(:);
new_ind = last_ind;
for iteration = 1:options.kmeans_max_iter %#ok<NASGU>
    distance = local_distance_to_centers(points, points(last_ind, :));
    [~, cluster_id] = min(distance, [], 2);
    for imu = 1:rank_mu
        cluster = find(cluster_id == imu);
        if isempty(cluster)
            new_ind(imu) = last_ind(imu);
            continue;
        end
        cluster_weight = weight(cluster);
        if sum(cluster_weight) <= 0
            cluster_weight = ones(size(cluster_weight));
        end
        centroid = sum(points(cluster, :) .* cluster_weight, 1) / ...
            sum(cluster_weight);
        distance_to_centroid = sum(abs( ...
            points(cluster, :) - centroid).^2, 2);
        [~, nearest] = min(distance_to_centroid);
        new_ind(imu) = cluster(nearest);
    end
    new_ind = local_unique_fill(new_ind, weight, rank_mu, ngrid);
    if isequal(new_ind(:), last_ind(:))
        break;
    end
    last_ind = new_ind(:);
end
ind_mu = reshape(new_ind, 1, []);
end

function points = local_grid_points(ngrid, options)
if isfield(options, 'points') && ~isempty(options.points)
    points = options.points;
    if size(points, 1) ~= ngrid
        error('ISDF:PointGridMismatch', ...
            'options.points must have one row per real-space grid point.');
    end
    return;
end
if isfield(options, 'fftgrid') && ~isempty(options.fftgrid)
    fftgrid = options.fftgrid(:).';
    if prod(fftgrid) ~= ngrid
        error('ISDF:GridSizeMismatch', ...
            ['prod(options.fftgrid) must match the number of ' ...
             'real-space grid points.']);
    end
    axes = cell(1, numel(fftgrid));
    for idimension = 1:numel(fftgrid)
        count = fftgrid(idimension);
        coordinate = (0:count-1) / count;
        coordinate(coordinate >= 0.5) = coordinate(coordinate >= 0.5) - 1;
        axes{idimension} = coordinate;
    end
    grids = cell(1, numel(fftgrid));
    [grids{:}] = ndgrid(axes{:});
    points = zeros(ngrid, numel(fftgrid));
    for idimension = 1:numel(fftgrid)
        points(:, idimension) = grids{idimension}(:);
    end
    return;
end
points = (0:ngrid-1).' / max(1, ngrid - 1);
end

% ---- K-means initialization and repair helpers ----

function ind_mu = local_kmeanspp_init(points, weight, rank_mu, options)
ngrid = size(points, 1);
ind_mu = zeros(rank_mu, 1);
ind_mu(1) = randi(ngrid);
for imu = 2:rank_mu
    centers = points(ind_mu(1:imu-1), :);
    distance = local_distance_to_centers(points, centers);
    min_distance = min(distance, [], 2);
    score = min_distance .* weight;
    score(ind_mu(1:imu-1)) = 0;
    if strcmpi(options.init, 'kmeans++')
        ind_mu(imu) = local_weighted_sample(score, 1, true);
    else
        [~, ind_mu(imu)] = max(score);
    end
end
end

function distance = local_distance_to_centers(points, centers)
distance = zeros(size(points, 1), size(centers, 1));
for icenter = 1:size(centers, 1)
    distance(:, icenter) = sum(abs(points - centers(icenter, :)).^2, 2);
end
end

function ind_mu = local_unique_fill(ind_mu, weight, rank_mu, ngrid)
ind_mu = ind_mu(:);
[~, first_position] = unique(ind_mu, 'stable');
duplicate_position = setdiff((1:numel(ind_mu)).', ...
    first_position, 'stable');
used = false(ngrid, 1);
used(ind_mu(first_position)) = true;
[~, order] = sort(weight(:), 'descend');
order = order(:);
for position = reshape(duplicate_position, 1, [])
    replacement = order(find(~used(order), 1, 'first'));
    ind_mu(position) = replacement;
    used(replacement) = true;
end
if numel(ind_mu) < rank_mu
    for position = (numel(ind_mu) + 1):rank_mu
        replacement = order(find(~used(order), 1, 'first'));
        ind_mu(position, 1) = replacement;
        used(replacement) = true;
    end
end
end

% ---- Weighted random sampling ----

function ind = local_weighted_sample(weight, nselect, replace)
weight = real(weight(:));
weight(~isfinite(weight) | weight < 0) = 0;
if all(weight == 0)
    weight = ones(size(weight));
end
ind = zeros(nselect, 1);
if replace
    total = sum(weight);
    for iselection = 1:nselect
        ind(iselection) = local_weighted_one(weight, total);
    end
else
    work_weight = weight;
    for iselection = 1:nselect
        total = sum(work_weight);
        if total <= 0
            remaining = find(work_weight >= 0);
            remaining = setdiff(remaining, ...
                ind(1:iselection-1), 'stable');
            ind(iselection) = remaining(1);
        else
            ind(iselection) = local_weighted_one(work_weight, total);
        end
        work_weight(ind(iselection)) = 0;
    end
end
end

function ind = local_weighted_one(weight, total)
cumulative = cumsum(weight);
target = rand * total;
ind = find(cumulative >= target, 1, 'first');
if isempty(ind)
    ind = numel(weight);
end
end
