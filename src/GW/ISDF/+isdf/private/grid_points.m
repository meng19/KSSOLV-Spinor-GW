function points = grid_points(ngrid, options)
%GRID_POINTS Coordinates used by weighted K-means point selection.

if isfield(options, 'points') && ~isempty(options.points)
    points = gather_if_gpu(options.points);
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
