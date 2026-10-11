function ind_mu = kmeanspp_init(points, weight, rank_mu, options)
%KMEANSPP_INIT Initialize K-means centers with weighted K-means++.

ngrid = size(points, 1);
ind_mu = zeros(rank_mu, 1);
ind_mu(1) = randi(ngrid);
for imu = 2:rank_mu
    centers = points(ind_mu(1:imu-1), :);
    [~, min_distance] = nearest_centers(points, centers, ...
        options.kmeans_distance_block_elements);
    score = min_distance .* weight;
    score(ind_mu(1:imu-1)) = 0;
    if strcmpi(options.init, 'kmeans++')
        ind_mu(imu) = weighted_sample(score, 1, true);
    else
        [~, ind_mu(imu)] = max(score);
    end
end
end
