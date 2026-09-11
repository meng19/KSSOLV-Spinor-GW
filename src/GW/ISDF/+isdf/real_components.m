function [components, cache_entries] = real_components( ...
    wfn, fft_template, idx, ispin, nspinor, bands, ...
    cache_entries, cache_all_bands, progress)
%ISDF.REAL_COMPONENTS Real-space spinor components with optional caching.

if nargin < 7 || isempty(cache_entries)
    cache_entries = cell(1, nspinor);
end
if nargin < 8
    cache_all_bands = false;
end
if nargin < 9
    progress = [];
end

bands = bands(:).';
components = cell(1, nspinor);
grid_size = size(fft_template);
for ispinor = 1:nspinor
    cached = cache_entries{ispinor};
    if local_cache_matches(cached, wfn, grid_size, bands)
        components{ispinor} = cached.values(:, bands);
        local_progress(progress, ispinor / nspinor, ...
            sprintf('real-space FFT reused spinor %d/%d', ispinor, nspinor));
        continue;
    end

    if cache_all_bands
        all_bands = 1:size(wfn.psi{ispin, ispinor}, 2);
        values = isdf.real_component( ...
            wfn, fft_template, idx, ispin, ispinor, all_bands, ...
            @(current, total) local_component_progress( ...
            progress, ispinor, nspinor, current, total));
        cache_entries{ispinor} = struct('values', values, ...
            'mill', wfn.mill, 'grid_size', grid_size);
        components{ispinor} = values(:, bands);
    else
        components{ispinor} = isdf.real_component( ...
            wfn, fft_template, idx, ispin, ispinor, bands, ...
            @(current, total) local_component_progress( ...
            progress, ispinor, nspinor, current, total));
    end
end
end

function local_component_progress(progress, ispinor, nspinor, current, total)
fraction = ((ispinor - 1) + current / total) / nspinor;
local_progress(progress, fraction, sprintf( ...
    'real-space FFT spinor %d/%d band %d/%d', ...
    ispinor, nspinor, current, total));
end

function local_progress(progress, fraction, stage)
if isa(progress, 'function_handle')
    progress(fraction, stage);
end
end

function tf = local_cache_matches(cached, wfn, grid_size, bands)
tf = isstruct(cached) && isfield(cached, 'values') && ...
    isequal(cached.grid_size, grid_size) && isequal(cached.mill, wfn.mill) && ...
    size(cached.values, 2) >= max(bands);
end
