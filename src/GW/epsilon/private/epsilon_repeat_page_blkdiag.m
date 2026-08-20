function repeated = epsilon_repeat_page_blkdiag(coeff, repeat_count)
%EPSILON_REPEAT_PAGE_BLKDIAG Repeat a reduced coefficient page block.

if repeat_count == 1
    repeated = coeff;
    return;
end
blocks = repmat({coeff}, 1, repeat_count);
repeated = epsilon_page_blkdiag(blocks);
end
