function epsilon_progress(block, local_work, message)
%EPSILON_PROGRESS Report global epsilon progress from inside a block.

if ~isfield(block, 'progress') || isempty(block.progress)
    return;
end
progress = block.progress;
local_work = min(max(local_work, 0), progress.block_work);
print_progress(progress.completed_before + local_work, ...
    progress.total_work, ...
    'Message', message, ...
    'Task', progress.task, ...
    'PercentStep', progress.percent_step);
end
