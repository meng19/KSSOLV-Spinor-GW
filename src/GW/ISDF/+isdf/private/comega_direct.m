function result = comega_direct(left, right, ev_occ, ev_unocc, freq, progress)
%COMEGA_DIRECT Direct frequency-page reduced polarizability.

if nargin < 6
    progress = [];
end
nmu = size(left{1}, 1);
nv = numel(ev_occ);
nc = numel(ev_unocc);
ncomponents = numel(left);
if ncomponents ~= numel(right)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end
result = complex(zeros(nmu, nmu, numel(freq), 'like', left{1}));
for iv = 1:nv
    for ic = 1:nc
        coefficient = complex(zeros(nmu, 1, 'like', left{1}));
        for icomponent = 1:ncomponents
            coefficient = coefficient + ...
                conj(left{icomponent}(:, iv)) .* right{icomponent}(:, ic);
        end
        energy_diff = ev_occ(iv) - ev_unocc(ic);
        for ifreq = 1:numel(freq)
            if freq(ifreq) == 0
                eden = 1 / energy_diff;
            else
                eden = energy_diff / (energy_diff^2 - freq(ifreq)^2);
            end
            result(:, :, ifreq) = result(:, :, ifreq) + ...
                coefficient * coefficient' * eden;
        end
        local_progress(progress, iv, ic, nv, nc);
    end
end
end

function local_progress(progress, iv, ic, nv, nc)
if isempty(progress)
    return;
end
pair_count = ic + (iv - 1) * nc;
npairs = max(1, nv * nc);
current = progress.completed_before + ...
    progress.block_work * pair_count / npairs;
left_band = local_band(progress, 'left_bands', iv);
right_band = local_band(progress, 'right_bands', ic);
print_progress(current, progress.total_work, ...
    'Message', sprintf('%s p%d/%d v%d c%d', ...
    progress.label, pair_count, npairs, left_band, right_band), ...
    'Task', progress.task, ...
    'PercentStep', progress.percent_step);
end

function band = local_band(progress, field, index)
if isfield(progress, field) && numel(progress.(field)) >= index
    band = progress.(field)(index);
else
    band = index;
end
end
