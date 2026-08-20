function ops = sigma_ops(ctx)
%SIGMA_OPS Select matrix construction and contraction handlers.

ops.name = ctx.method;
switch ctx.method
    case 'direct'
        ops.matrix_elements = @(block) ...
            sigma_matrix_elements(ctx, block, false);
        ops.contract = @(block, matrix_elements) ...
            sigma_contract_full(ctx, block, matrix_elements);
    case 'matrix_elements'
        ops.matrix_elements = @(block) ...
            sigma_matrix_elements(ctx, block, true);
        ops.contract = @(block, matrix_elements) ...
            sigma_contract_full(ctx, block, matrix_elements);
    case 'reduced_basis'
        ops.matrix_elements = @(block) ...
            sigma_matrix_elements(ctx, block, true);
        ops.contract = @(block, matrix_elements) ...
            sigma_contract_reduced(ctx, block, matrix_elements);
end
end
