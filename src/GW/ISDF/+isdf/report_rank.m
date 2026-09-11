function report_rank(options, ngrid, nleft, nright)
%ISDF.REPORT_RANK Print each distinct ISDF-space rank once per calculation.
%   ISDF.REPORT_RANK('reset') starts a new reporting scope.

persistent printed_keys;
if ischar(options) || isstring(options)
    if strcmpi(char(options), 'reset')
        printed_keys = {};
        return;
    end
    error('ISDF:RankReportAction', 'The only rank-report action is ''reset''.');
end
if isfield(options, 'print_rank') && ~options.print_rank
    return;
end
if isempty(printed_keys)
    printed_keys = {};
end

product_type = 'generic';
if isfield(options, 'product_type') && ~isempty(options.product_type)
    product_type = upper(char(options.product_type));
end
key = sprintf('%s:%s:%s:%d:%d:%d:%d:%d:%d:%.16g', ...
    product_type, lower(char(options.sample_method)), ...
    lower(char(options.rank_source)), ngrid, nleft, nright, options.rank, ...
    options.recommended_rank, options.max_rank, options.rank_ratio);
if any(strcmp(printed_keys, key))
    return;
end
printed_keys{end + 1} = key;

fprintf(['\nISDF %s space: %d interpolation points for %d left x %d right ' ...
    'bands; recommended ceil(sqrt(%d*%d)*%.3g) = %d\n'], ...
    product_type, options.rank, nleft, nright, nleft, nright, ...
    options.rank_ratio, options.recommended_rank);
end
