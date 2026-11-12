function [eps, sys, options, syms] = ...
    sigma_checkpoint_inputs(eps, sig, sys, options, syms)
%SIGMA_CHECKPOINT_INPUTS Resolve optional epsilon checkpoint inputs.

checkpoint_file = local_checkpoint_file(eps, sig);
if ~isempty(checkpoint_file)
    [eps, checkpoint_sys, checkpoint_options, checkpoint_syms] = ...
        epsilon_checkpoint_load(checkpoint_file);
    if isempty(sys)
        sys = checkpoint_sys;
    end
    if isempty(options)
        options = checkpoint_options;
    end
    if isempty(syms)
        syms = checkpoint_syms;
    end
end

if isempty(sys) || isempty(options) || isempty(syms)
    error('GW:SigmaCheckpointInputs', ...
        ['Provide sys, options, and syms, or pass an epsilon checkpoint ' ...
         'that was written by epsilon.']);
end
end

function checkpoint_file = local_checkpoint_file(eps, sig)
checkpoint_file = '';
if ischar(eps) || (isstring(eps) && isscalar(eps))
    checkpoint_file = eps;
elseif isempty(eps) && isstruct(sig) && ...
        isfield(sig, 'epsilon_checkpoint_file') && ...
        ~isempty(sig.epsilon_checkpoint_file)
    checkpoint_file = sig.epsilon_checkpoint_file;
elseif isempty(eps) && isstruct(sig) && ...
        isfield(sig, 'epsilon_checkpoint') && ...
        ~isempty(sig.epsilon_checkpoint)
    checkpoint_file = sig.epsilon_checkpoint;
end
end
