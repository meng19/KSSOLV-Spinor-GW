function work = gw_block_work(block, fallback_work)
%GW_BLOCK_WORK Return the assigned progress work, or a caller-provided fallback.

if isfield(block, 'progress') && isfield(block.progress, 'block_work') && ...
        ~isempty(block.progress.block_work)
    work = block.progress.block_work;
else
    work = fallback_work;
end
end
