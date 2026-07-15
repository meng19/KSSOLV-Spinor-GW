function ind_mu = isdf_kmeans_indices(weight, options)
%ISDF_KMEANS_INDICES Weighted K-means interpolation point selection.

rank_mu = options.rank;
weight = real(weight(:));
points = isdf_grid_points(numel(weight), options);
[ngrid, ~] = size(points);

if rank_mu > ngrid
    error('ISDF:RankTooLarge', 'ISDF rank cannot exceed the number of grid points.');
end

switch lower(options.init)
    case 'sequence'
        ind_mu = round(linspace(1, ngrid, rank_mu)).';
    case 'random'
        ind_mu = randperm(ngrid, rank_mu).';
    case 'wrs'
        ind_mu = isdf_weighted_sample(abs(weight), rank_mu, false);
    case {'kmeans++', 'kmeans++_false'}
        ind_mu = local_kmeanspp_init(points, abs(weight), rank_mu, options);
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
for iter = 1:options.kmeans_max_iter
    dist = local_distance_to_centers(points, points(last_ind, :));
    [~, cluster_id] = min(dist, [], 2);

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

        centroid = sum(points(cluster, :) .* cluster_weight, 1) / sum(cluster_weight);
        dist_to_centroid = sum(abs(points(cluster, :) - centroid).^2, 2);
        [~, nearest] = min(dist_to_centroid);
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

function ind_mu = local_kmeanspp_init(points, weight, rank_mu, options)
ngrid = size(points, 1);
ind_mu = zeros(rank_mu, 1);
ind_mu(1) = randi(ngrid);

for imu = 2:rank_mu
    centers = points(ind_mu(1:imu-1), :);
    dist = local_distance_to_centers(points, centers);
    min_dist = min(dist, [], 2);
    score = min_dist .* weight;
    score(ind_mu(1:imu-1)) = 0;
    if strcmpi(options.init, 'kmeans++')
        ind_mu(imu) = isdf_weighted_sample(score, 1, true);
    else
        [~, ind_mu(imu)] = max(score);
    end
end
end

function dist = local_distance_to_centers(points, centers)
dist = zeros(size(points, 1), size(centers, 1));
for icenter = 1:size(centers, 1)
    dist(:, icenter) = sum(abs(points - centers(icenter, :)).^2, 2);
end
end

function ind_mu = local_unique_fill(ind_mu, weight, rank_mu, ngrid)
ind_mu = ind_mu(:);
[~, first_pos] = unique(ind_mu, 'stable');
duplicate_pos = setdiff((1:numel(ind_mu)).', first_pos, 'stable');

used = false(ngrid, 1);
used(ind_mu(first_pos)) = true;
[~, order] = sort(weight(:), 'descend');
order = order(:);

for ipos = reshape(duplicate_pos, 1, [])
    replacement = order(find(~used(order), 1, 'first'));
    ind_mu(ipos) = replacement;
    used(replacement) = true;
end

if numel(ind_mu) < rank_mu
    for ipos = (numel(ind_mu) + 1):rank_mu
        replacement = order(find(~used(order), 1, 'first'));
        ind_mu(ipos, 1) = replacement;
        used(replacement) = true;
    end
end
end
