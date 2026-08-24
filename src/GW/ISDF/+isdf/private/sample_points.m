function [ind_mu, products] = sample_points(left, right, options)
%SAMPLE_POINTS Select ISDF interpolation points for component products.

if ~iscell(left)
    left = {left};
end
if ~iscell(right)
    right = {right};
end
rng(options.seed, 'twister');
products = [];

switch lower(options.sample_method)
    case 'qrcp'
        products = component_products(left, right, [], []);
        ind_mu = qrcp_sample(products, options.rank);
    case {'qrcp_randomized', 'randomized_qrcp', 'default'}
        if numel(left) == 1
            ind_mu = scalar_randomized_sample( ...
                left{1}, right{1}, options);
        else
            npairs = size(left{1}, 2) * size(right{1}, 2);
            projection_rank = max(options.rank, ...
                ceil(options.random_oversampling * options.rank));
            projection_rank = min(projection_rank, npairs);
            projection = randn_like(npairs, projection_rank, left{1});
            has_complex = any(cellfun(@(x) ~isreal(x), left)) || ...
                any(cellfun(@(x) ~isreal(x), right));
            if has_complex
                projection = projection + 1i * ...
                    randn_like(npairs, projection_rank, left{1});
            end
            compressed_products = component_products( ...
                left, right, [], projection);
            ind_mu = qrcp_sample(compressed_products, options.rank);
        end
    case 'kmeans'
        if numel(left) == 1
            weight = scalar_weight(left{1}, right{1}, options);
        else
            weight = component_weight(left, right, options);
        end
        ind_mu = kmeans_sample(weight, options);
    otherwise
        error('ISDF:UnknownSampleMethod', ...
            ['Unknown ISDF sample_method "%s". Supported methods: qrcp, ' ...
             'qrcp_randomized, kmeans, default.'], options.sample_method);
end
end
