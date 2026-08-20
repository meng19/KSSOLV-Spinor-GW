function ops = epsilon_ops(ctx)
%EPSILON_OPS Select epsilon handlers by workflow.

ops.name = ctx.method;
switch ctx.method
    case 'direct'
        ops.init = @(iq) epsilon_full_init(ctx, iq);
        ops.evaluate = @(block) epsilon_full_evaluate(ctx, block, false);
        ops.accumulate = @(acc, contribution, block) ...
            epsilon_full_accumulate(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) ...
            epsilon_full_finalize(ctx, eps, acc, iq);
    case 'matrix_elements'
        ops.init = @(iq) epsilon_full_init(ctx, iq);
        ops.evaluate = @(block) epsilon_full_evaluate(ctx, block, true);
        ops.accumulate = @(acc, contribution, block) ...
            epsilon_full_accumulate(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) ...
            epsilon_full_finalize(ctx, eps, acc, iq);
    case 'reduced_basis'
        output_mode = lower(ctx.eps.isdf.output);
        need_full_inverse = any(strcmp( ...
            output_mode, {'full_inverse', 'both'}));
        need_screened_w = any(strcmp( ...
            output_mode, {'screened_w', 'both'}));
        if ~need_full_inverse && ~need_screened_w
            error('ISDF:ReducedEpsilonOutput', ...
                'Unknown ISDF reduced-basis epsilon output "%s".', ...
                ctx.eps.isdf.output);
        end
        ops.init = @(iq) epsilon_reduced_init( ...
            ctx, iq, need_full_inverse, need_screened_w);
        ops.evaluate = @(block) epsilon_reduced_evaluate(ctx, block);
        ops.accumulate = @(acc, contribution, block) ...
            epsilon_reduced_accumulate(ctx, acc, contribution, block);
        ops.finalize = @(eps, acc, iq) ...
            epsilon_reduced_finalize(ctx, eps, acc, iq);
end
end
