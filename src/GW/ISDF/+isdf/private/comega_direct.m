function result = comega_direct(left, right, ev_occ, ev_unocc, freq)
%COMEGA_DIRECT Direct frequency-page reduced polarizability.

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
    end
end
end
