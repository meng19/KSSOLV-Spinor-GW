function ind_mu = kmeans_sample(weight, options)
%KMEANS_SAMPLE Select interpolation points using weighted K-means.

rank_mu = options.rank;
weight = gather_if_gpu(real(weight(:)));
points = grid_points(numel(weight), options);
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
        ind_mu = weighted_sample(abs(weight), rank_mu, false);
    case {'kmeans++', 'kmeans++_false'}
        ind_mu = kmeanspp_init(points, abs(weight), rank_mu, options);
    otherwise
        error('ISDF:UnknownKmeansInit', ...
            'Unknown ISDF K-means init "%s".', options.init);
end

ind_mu = unique_fill(ind_mu, weight, rank_mu, ngrid);
if options.kmeans_max_iter == 0
    ind_mu = reshape(ind_mu, 1, []);
    return;
end

last_ind = ind_mu(:);
new_ind = last_ind;
for iteration = 1:options.kmeans_max_iter %#ok<NASGU>
    [cluster_id, ~] = nearest_centers(points, points(last_ind, :), ...
        options.kmeans_distance_block_elements);
    % Group points by cluster once instead of running one O(npoints)
    % find() per center; large grids spend most of the time there.
    [sorted_id, order] = sort(cluster_id);
    counts = accumarray(sorted_id, 1, [rank_mu, 1]);
    first_pos = cumsum([1; counts(1:end - 1)]);
    for imu = 1:rank_mu
        cluster = order(first_pos(imu):first_pos(imu) + counts(imu) - 1);
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
    new_ind = unique_fill(new_ind, weight, rank_mu, ngrid);
    if isequal(new_ind(:), last_ind(:))
        break;
    end
    last_ind = new_ind(:);
end
ind_mu = reshape(new_ind, 1, []);
end
