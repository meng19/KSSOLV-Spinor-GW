function gw_block_progress(block, local_work, message)
%GW_BLOCK_PROGRESS Report global progress from a prepared epsilon/sigma block.

if ~isfield(block, 'progress') || isempty(block.progress)
    return;
end
progress = block.progress;
local_work = min(max(local_work, 0), progress.block_work);
update_interval = 5;
if isfield(progress, 'update_interval') && ~isempty(progress.update_interval)
    update_interval = progress.update_interval;
end
print_progress(progress.completed_before + local_work, progress.total_work, ...
    'Message', message, 'Task', progress.task, ...
    'PercentStep', progress.percent_step, 'UpdateInterval', update_interval);
end
