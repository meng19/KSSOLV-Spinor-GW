function weight = component_weight(left, right, options)
%COMPONENT_WEIGHT Compute component-product sampling weights.

if nargin < 3 || isempty(options)
    options = struct();
end
if ~isfield(options, 'weight_was_set') || ~options.weight_was_set
    options.weight = 'prod';
elseif ~isfield(options, 'weight') || isempty(options.weight)
    options.weight = 'prod';
end

ngrid = size(left{1}, 1);
nleft = size(left{1}, 2);
nright = size(right{1}, 2);
local_validate_components(left, right, ngrid, nleft, nright);

switch lower(options.weight)
    case {'prod', 'multiply'}
        weight = local_product_weight(left, right, ngrid);
    case 'add'
        weight = local_density_weight(left, right, ngrid);
    case 'power'
        if ~isfield(options, 'power') || isempty(options.power)
            options.power = 1;
        end
        weight = local_product_weight(left, right, ngrid) .^ ...
            (options.power / 2);
    otherwise
        error('ISDF:UnknownWeight', ...
            ['Unknown ISDF weight "%s". Supported weights: ' ...
             'prod, add, power.'], options.weight);
end
weight = gather_if_gpu(real(weight(:)));
weight(~isfinite(weight) | weight < 0) = 0;
if all(weight == 0)
    weight = ones(size(weight));
end
end

function local_validate_components(left, right, ngrid, nleft, nright)
for icomponent = 1:numel(left)
    left_i = left{icomponent};
    right_i = right{icomponent};
    if size(left_i, 1) ~= ngrid || size(right_i, 1) ~= ngrid || ...
            size(left_i, 2) ~= nleft || size(right_i, 2) ~= nright
        error('ISDF:ComponentSizeMismatch', ...
            'All components must share grid and band dimensions.');
    end
end
end

function weight = local_product_weight(left, right, ngrid)
weight = zeros(ngrid, 1, 'like', left{1});
for icomponent = 1:numel(left)
    left_i = left{icomponent};
    right_i = right{icomponent};
    for jcomponent = 1:numel(left)
        left_j = left{jcomponent};
        right_j = right{jcomponent};
        weight = weight + sum(conj(left_i) .* left_j, 2) .* ...
            sum(right_i .* conj(right_j), 2);
    end
end
end

function weight = local_density_weight(left, right, ngrid)
weight = zeros(ngrid, 1, 'like', left{1});
for icomponent = 1:numel(left)
    weight = weight + sum(abs(left{icomponent}).^2, 2) + ...
        sum(abs(right{icomponent}).^2, 2);
end
end
