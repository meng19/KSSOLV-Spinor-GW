function distance = distance_to_centers(points, centers)
%DISTANCE_TO_CENTERS Squared distance from points to each center.

distance = zeros(size(points, 1), size(centers, 1));
for icenter = 1:size(centers, 1)
    distance(:, icenter) = sum(abs(points - centers(icenter, :)).^2, 2);
end
end
