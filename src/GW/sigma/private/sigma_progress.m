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
    'PercentStep', progress.percent_step);
end
