function ind = weighted_one(weight, total)
%WEIGHTED_ONE Draw one weighted index.

cumulative = cumsum(weight);
target = rand * total;
ind = find(cumulative >= target, 1, 'first');
if isempty(ind)
    ind = numel(weight);
end
end
