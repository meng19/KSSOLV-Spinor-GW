function epsilon_checkpoint_save(filename, eps, request_eps, sys, options, syms)
%EPSILON_CHECKPOINT_SAVE Persist epsilon output for a later sigma run.

filename = char(filename);
[folder, ~, ~] = fileparts(filename);
if ~isempty(folder) && exist(folder, 'dir') ~= 7
    error('GW:EpsilonCheckpointDirectory', ...
        'Checkpoint directory does not exist: %s', folder);
end

checkpoint = struct();
checkpoint.version = 1;
checkpoint.eps = eps;
checkpoint.request_eps = request_eps;
checkpoint.sys = sys;
checkpoint.options = options;
checkpoint.syms = syms;
save(filename, 'checkpoint', '-v7.3');
fprintf('Saved epsilon checkpoint: %s\n', filename);
end
