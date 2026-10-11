function [cluster_id, min_distance] = nearest_centers( ...
        points, centers, block_elements)
%NEAREST_CENTERS Assign points without materializing all point-center distances.

if nargin < 3 || isempty(block_elements)
    block_elements = 2e7;
end
npoints = size(points, 1);
ncenters = size(centers, 1);
if ncenters < 1
    error('ISDF:EmptyKmeansCenters', 'At least one K-means center is required.');
end

% DISTANCE_TO_CENTERS creates block_rows*ncenters real values.  Bound that
% temporary explicitly; the output vectors have only npoints entries.
block_rows = max(1, floor(block_elements / ncenters));
cluster_id = zeros(npoints, 1);
min_distance = zeros(npoints, 1, 'like', real(points(:, 1)));
for first_point = 1:block_rows:npoints
    last_point = min(npoints, first_point + block_rows - 1);
    block = first_point:last_point;
    distance = distance_to_centers(points(block, :), centers);
    [min_distance(block), cluster_id(block)] = min(distance, [], 2);
end
end
