function distance = distance_to_centers(points, centers)
%DISTANCE_TO_CENTERS Squared distance from points to each center.
%
% Use ||x-c||^2 = ||x||^2 + ||c||^2 - 2 Re(x*c') so that the
% point-center distance matrix is formed by one BLAS matrix multiply.
% This avoids allocating one ngrid-by-ndimension temporary per center in
% every K-means iteration.

point_norm = sum(abs(points).^2, 2);
center_norm = sum(abs(centers).^2, 2).';
distance = point_norm + center_norm - 2 * real(points * centers');

% Round-off can make distances between identical points slightly negative.
distance = max(distance, 0);
end
