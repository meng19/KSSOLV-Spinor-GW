function acc = epsilon_full_accumulate(ctx, acc, contribution, block)
%EPSILON_FULL_ACCUMULATE Accumulate full matrix-element contributions.

if isfield(contribution, 'stream_direct')
    acc = epsilon_direct_stream_accumulate(ctx, acc, block);
    return;
end

if ctx.save_mem
    acc.chi0 = epsilon_add_state_batch( ...
        acc.chi0, contribution.gme, contribution.eden, block.g_maps);
    return;
end

gme = epsilon_mapped_gme(contribution.gme, block.g_maps);
eden = epsilon_repeat_eden(contribution.eden, numel(block.g_maps));
if isempty(acc.deferred{block.ispin, block.ik_fbz})
    acc.deferred{block.ispin, block.ik_fbz} = {{gme, eden}};
else
    acc.deferred{block.ispin, block.ik_fbz}{end + 1} = {gme, eden};
end
end
