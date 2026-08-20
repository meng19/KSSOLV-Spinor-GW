function inverse_pages = epsilon_invert_pages(nmtx, coulg, chi0)
%EPSILON_INVERT_PAGES Invert each dynamic epsilon page independently.

if isa(chi0, 'gpuArray')
    if ~isa(coulg, 'gpuArray')
        coulg = gpuArray(coulg);
    end
    identity = eye(nmtx, 'like', chi0);
    epsilon_pages = repmat(identity, 1, 1, size(chi0, 3));
    epsilon_pages = epsilon_pages - bsxfun(@times, coulg(:), chi0);
    inverse_pages = gather(pagefun(@inv, epsilon_pages));
    return;
end
identity = repmat(eye(nmtx), 1, 1, size(chi0, 3));
epsilon_pages = identity - bsxfun(@times, coulg(:), chi0);
inverse_cells = cell(1, size(epsilon_pages, 3));
for ifreq = 1:size(epsilon_pages, 3)
    inverse_cells{ifreq} = inv(epsilon_pages(:, :, ifreq));
end
inverse_pages = cat(3, inverse_cells{:});
end
