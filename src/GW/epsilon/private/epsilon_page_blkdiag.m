function combined = epsilon_page_blkdiag(blocks)
%EPSILON_PAGE_BLKDIAG Page-wise block diagonal assembly.

if numel(blocks) == 1
    combined = blocks{1};
    return;
end

nfreq = size(blocks{1}, 3);
total_rows = 0;
total_cols = 0;
for iblock = 1:numel(blocks)
    total_rows = total_rows + size(blocks{iblock}, 1);
    total_cols = total_cols + size(blocks{iblock}, 2);
end
% Fill the output pages directly.  Holding every completed page in a cell
% array and then calling cat(3, ...) temporarily keeps two full copies.
combined = complex(zeros(total_rows, total_cols, nfreq, ...
    'like', blocks{1}));
for ifreq = 1:nfreq
    row_start = 1;
    col_start = 1;
    for iblock = 1:numel(blocks)
        block = blocks{iblock}(:, :, ifreq);
        row_idx = row_start:(row_start + size(block, 1) - 1);
        col_idx = col_start:(col_start + size(block, 2) - 1);
        combined(row_idx, col_idx, ifreq) = block;
        row_start = row_idx(end) + 1;
        col_start = col_idx(end) + 1;
    end
end
end
