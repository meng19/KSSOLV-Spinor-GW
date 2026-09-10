function sigma_progress(block, local_work, message)
%SIGMA_PROGRESS Report global sigma progress from inside a q block.

if ~isfield(block, 'progress') || isempty(block.progress)
    return;
end
progress = block.progress;
local_work = min(max(local_work, 0), progress.block_work);
print_progress(progress.completed_before + local_work, ...
    progress.total_work, ...
    'Message', message, ...
    'Task', progress.task, ...
    'PercentStep', progress.percent_step, ...
    'UpdateInterval', local_update_interval(progress));
end

function interval = local_update_interval(progress)
interval = 5;
if isfield(progress, 'update_interval') && ~isempty(progress.update_interval)
    interval = progress.update_interval;
end
end
