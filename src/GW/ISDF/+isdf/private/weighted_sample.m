function ind = weighted_sample(weight, nselect, replace)
%WEIGHTED_SAMPLE Draw weighted sample indices.

weight = gather_if_gpu(real(weight(:)));
weight(~isfinite(weight) | weight < 0) = 0;
if all(weight == 0)
    weight = ones(size(weight));
end
ind = zeros(nselect, 1);
if replace
    total = sum(weight);
    for iselection = 1:nselect
        ind(iselection) = weighted_one(weight, total);
    end
else
    work_weight = weight;
    for iselection = 1:nselect
        total = sum(work_weight);
        if total <= 0
            remaining = find(work_weight >= 0);
            remaining = setdiff(remaining, ...
                ind(1:iselection-1), 'stable');
            ind(iselection) = remaining(1);
        else
            ind(iselection) = weighted_one(work_weight, total);
        end
        work_weight(ind(iselection)) = 0;
    end
end
end
