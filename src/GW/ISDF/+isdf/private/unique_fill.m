function ind_mu = unique_fill(ind_mu, weight, rank_mu, ngrid)
%UNIQUE_FILL Replace duplicate or missing sample indices by weighted order.

ind_mu = ind_mu(:);
[~, first_position] = unique(ind_mu, 'stable');
duplicate_position = setdiff((1:numel(ind_mu)).', ...
    first_position, 'stable');
used = false(ngrid, 1);
used(ind_mu(first_position)) = true;
[~, order] = sort(weight(:), 'descend');
order = order(:);
for position = reshape(duplicate_position, 1, [])
    replacement = order(find(~used(order), 1, 'first'));
    ind_mu(position) = replacement;
    used(replacement) = true;
end
if numel(ind_mu) < rank_mu
    for position = (numel(ind_mu) + 1):rank_mu
        replacement = order(find(~used(order), 1, 'first'));
        ind_mu(position, 1) = replacement;
        used(replacement) = true;
    end
end
end
