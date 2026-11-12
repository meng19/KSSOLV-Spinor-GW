function [eps, sys, options, syms, request_eps] = epsilon_checkpoint_load(filename)
%EPSILON_CHECKPOINT_LOAD Load epsilon output and its calculation inputs.

filename = char(filename);
if exist(filename, 'file') ~= 2
    error('GW:EpsilonCheckpointMissing', ...
        'Epsilon checkpoint does not exist: %s', filename);
end

loaded = load(filename, 'checkpoint');
if ~isfield(loaded, 'checkpoint') || ~isstruct(loaded.checkpoint)
    error('GW:EpsilonCheckpointFormat', ...
        'File is not an epsilon checkpoint: %s', filename);
end
checkpoint = loaded.checkpoint;
required = {'eps', 'sys', 'options', 'syms'};
if ~all(isfield(checkpoint, required))
    error('GW:EpsilonCheckpointFormat', ...
        'Epsilon checkpoint is missing required data: %s', filename);
end

eps = checkpoint.eps;
sys = checkpoint.sys;
options = checkpoint.options;
syms = checkpoint.syms;
if isfield(checkpoint, 'request_eps')
    request_eps = checkpoint.request_eps;
else
    request_eps = [];
end
fprintf('Loaded epsilon checkpoint: %s\n', filename);
end
