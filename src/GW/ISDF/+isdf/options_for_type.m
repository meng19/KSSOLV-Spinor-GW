function options = options_for_type(options, product_type)
%ISDF.OPTIONS_FOR_TYPE Apply optional product-space-specific rank settings.
%
% Supported product_type values are 'vc', 'vn', and 'nn'.  A specific
% rank_<type> takes precedence over rank; a rank_ratio_<type> takes
% precedence over rank_ratio and requests automatic rank selection from
% the corresponding product-space dimensions.

product_type = lower(char(product_type));
if ~any(strcmp(product_type, {'vc', 'vn', 'nn'}))
    error('ISDF:UnknownProductType', ...
        'Unknown ISDF product-space type "%s".', product_type);
end

rank_field = ['rank_' product_type];
ratio_field = ['rank_ratio_' product_type];
if isfield(options, rank_field) && ~isempty(options.(rank_field))
    options.rank = options.(rank_field);
elseif isfield(options, ratio_field) && ~isempty(options.(ratio_field))
    options.rank_ratio = options.(ratio_field);
    if isfield(options, 'rank')
        options = rmfield(options, 'rank');
    end
end
end
