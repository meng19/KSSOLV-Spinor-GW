function combined = epsilon_page_blkdiag(blocks)
%EPSILON_PAGE_BLKDIAG Page-wise block diagonal assembly.

nfreq = size(blocks{1}, 3);
page_cells = cell(1, nfreq);
for ifreq = 1:nfreq
    total_rows = 0;
    total_cols = 0;
    for iblock = 1:numel(blocks)
        total_rows = total_rows + size(blocks{iblock}, 1);
        total_cols = total_cols + size(blocks{iblock}, 2);
    end
    page = complex(zeros(total_rows, total_cols, 'like', blocks{1}));
    row_start = 1;
    col_start = 1;
    for iblock = 1:numel(blocks)
        block = blocks{iblock}(:, :, ifreq);
        row_idx = row_start:(row_start + size(block, 1) - 1);
        col_idx = col_start:(col_start + size(block, 2) - 1);
        page(row_idx, col_idx) = block;
        row_start = row_idx(end) + 1;
        col_start = col_idx(end) + 1;
    end
    page_cells{ifreq} = page;
end
combined = cat(3, page_cells{:});
end
