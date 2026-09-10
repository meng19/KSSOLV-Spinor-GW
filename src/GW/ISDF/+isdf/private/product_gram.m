function [c1, c2] = product_gram(left_components, right_components, ind_mu)
%PRODUCT_GRAM Build PP_mu' and P_mu P_mu' without forming all products.

if ~iscell(left_components)
    left_components = {left_components};
end
if ~iscell(right_components)
    right_components = {right_components};
end
if numel(left_components) ~= numel(right_components)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end

ngrid = size(left_components{1}, 1);
nleft = size(left_components{1}, 2);
nright = size(right_components{1}, 2);
for icomponent = 1:numel(left_components)
    if size(left_components{icomponent}, 1) ~= ngrid || ...
            size(right_components{icomponent}, 1) ~= ngrid
        error('ISDF:GridMismatch', ...
            'All left and right components must share the same grid size.');
    end
    if size(left_components{icomponent}, 2) ~= nleft || ...
            size(right_components{icomponent}, 2) ~= nright
        error('ISDF:BandMismatch', ...
            'All left and right components must share the same band counts.');
    end
end

ind_mu = ind_mu(:);
nmu = numel(ind_mu);
c1 = complex(zeros(ngrid, nmu, 'like', left_components{1}));
c2 = complex(zeros(nmu, nmu, 'like', left_components{1}));
left_mu = cell(size(left_components));
right_mu = cell(size(right_components));
for icomponent = 1:numel(left_components)
    left_mu{icomponent} = left_components{icomponent}(ind_mu, :);
    right_mu{icomponent} = right_components{icomponent}(ind_mu, :);
end
for icomponent = 1:numel(left_components)
    left_i = left_components{icomponent};
    right_i = right_components{icomponent};
    left_i_mu = left_mu{icomponent};
    right_i_mu = right_mu{icomponent};
    for jcomponent = 1:numel(left_components)
        left_j_mu = left_mu{jcomponent};
        right_j_mu = right_mu{jcomponent};
        c1 = c1 + (conj(left_i) * left_j_mu.') .* ...
            (right_i * conj(right_j_mu).');
        c2 = c2 + (conj(left_i_mu) * left_j_mu.') .* ...
            (right_i_mu * conj(right_j_mu).');
    end
end
end
