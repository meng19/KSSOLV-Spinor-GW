function points = isdf_grid_points(ngrid, options)
%ISDF_GRID_POINTS Return real-space coordinates used by K-means sampling.

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
            'prod(options.fftgrid) must match the number of real-space grid points.');
    end

    axes = cell(1, numel(fftgrid));
    for idim = 1:numel(fftgrid)
        n = fftgrid(idim);
        x = (0:n-1) / n;
        x(x >= 0.5) = x(x >= 0.5) - 1;
        axes{idim} = x;
    end

    grids = cell(1, numel(fftgrid));
    [grids{:}] = ndgrid(axes{:});
    points = zeros(ngrid, numel(fftgrid));
    for idim = 1:numel(fftgrid)
        points(:, idim) = grids{idim}(:);
    end
    return;
end

points = ((0:ngrid-1).' / max(1, ngrid - 1));
end
