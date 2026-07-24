function space = build_space(left, right, idx_q, fftgrid, options)
if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
space = isdf_build_space(left, right, idx_q, fftgrid, options);
end
