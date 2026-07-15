function ind = isdf_weighted_sample(weight, nselect, replace)
%ISDF_WEIGHTED_SAMPLE Draw weighted indices without requiring toolboxes.

if nargin < 3
    replace = false;
end

weight = real(weight(:));
weight(~isfinite(weight) | weight < 0) = 0;
if all(weight == 0)
    weight = ones(size(weight));
end

if replace
    ind = zeros(nselect, 1);
    total = sum(weight);
    for ii = 1:nselect
        ind(ii) = local_weighted_one(weight, total);
    end
else
    ind = zeros(nselect, 1);
    work_weight = weight;
    for ii = 1:nselect
        total = sum(work_weight);
        if total <= 0
            remaining = find(work_weight >= 0);
            remaining = setdiff(remaining, ind(1:ii-1), 'stable');
            ind(ii) = remaining(1);
        else
            ind(ii) = local_weighted_one(work_weight, total);
        end
        work_weight(ind(ii)) = 0;
    end
end
end

function ind = local_weighted_one(weight, total)
cumulative = cumsum(weight);
target = rand * total;
ind = find(cumulative >= target, 1, 'first');
if isempty(ind)
    ind = numel(weight);
end
end
