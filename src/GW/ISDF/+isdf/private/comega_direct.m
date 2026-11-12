function result = comega_direct(left, right, ev_occ, ev_unocc, ...
    freq, progress, options)
%COMEGA_DIRECT Exact frequency pages using blocked pair-matrix GEMMs.
% Each P chunk is constructed once and reused for every frequency page.

if nargin < 6
    progress = [];
end
if nargin < 7
    options = struct();
end
if ~isfield(options, 'cauchy_coefficient_block_bytes')
    options.cauchy_coefficient_block_bytes = 512 * 1024^2;
end
nmu = size(left{1}, 1);
nv = numel(ev_occ);
nc = numel(ev_unocc);
ncomponents = numel(left);
if ncomponents ~= numel(right)
    error('ISDF:ComponentMismatch', ...
        'Left and right component counts must match.');
end
chunks = local_coefficient_chunks(left, right, options);
result = complex(zeros(nmu, nmu, numel(freq), 'like', left{1}));
ev_occ_work = ev_occ;
ev_unocc_work = ev_unocc;
if isa(left{1}, 'gpuArray')
    ev_occ_work = gpuArray(ev_occ_work);
    ev_unocc_work = gpuArray(ev_unocc_work);
end
for ichunk = 1:numel(chunks)
    chunk = chunks(ichunk);
    P = local_chunk_coefficients(chunk, left, right);
    columns = chunk.first:chunk.last;
    iv = mod(columns - 1, nv) + 1;
    ic = floor((columns - 1) / nv) + 1;
    energy_diff = reshape( ...
        ev_occ_work(iv) - ev_unocc_work(ic), 1, []);
    for ifreq = 1:numel(freq)
        if freq(ifreq) == 0
            eden = 1 ./ energy_diff;
        else
            eden = energy_diff ./ (energy_diff.^2 - freq(ifreq)^2);
        end
        result(:, :, ifreq) = result(:, :, ifreq) + ...
            (P .* eden) * P';
    end
    local_progress(progress, chunk.last, nv, nc);
end
end

function chunks = local_coefficient_chunks(left, right, options)
% Describe pair-column chunks.  Each is used once for all frequency pages,
% so caching every P chunk would only increase memory without saving work.
nmu = size(left{1}, 1);
npairs = size(left{1}, 2) * size(right{1}, 2);
pair_block = max(1, floor(options.cauchy_coefficient_block_bytes / ...
    max(1, 32 * nmu)));
pair_block = min(npairs, pair_block);
chunks = struct('first', {}, 'last', {});
for first = 1:pair_block:npairs
    last = min(npairs, first + pair_block - 1);
    chunks(end + 1) = struct('first', first, 'last', last); %#ok<AGROW>
end
end

function P = local_chunk_coefficients(chunk, left, right)
nmu = size(left{1}, 1);
nv = size(left{1}, 2);
columns = chunk.first:chunk.last;
iv = mod(columns - 1, nv) + 1;
ic = floor((columns - 1) / nv) + 1;
P = complex(zeros(nmu, numel(columns), 'like', left{1}));
for icomponent = 1:numel(left)
    P = P + conj(left{icomponent}(:, iv)) .* right{icomponent}(:, ic);
end
end

function local_progress(progress, pair_count, nv, nc)
if isempty(progress)
    return;
end
iv = mod(pair_count - 1, nv) + 1;
ic = floor((pair_count - 1) / nv) + 1;
npairs = max(1, nv * nc);
current = progress.completed_before + ...
    progress.block_work * pair_count / npairs;
left_band = local_band(progress, 'left_bands', iv);
right_band = local_band(progress, 'right_bands', ic);
print_progress(current, progress.total_work, ...
    'Message', sprintf('%s p%d/%d v%d c%d', ...
    progress.label, pair_count, npairs, left_band, right_band), ...
    'Task', progress.task, ...
    'PercentStep', progress.percent_step, ...
    'UpdateInterval', local_update_interval(progress));
end

function interval = local_update_interval(progress)
interval = 5;
if isfield(progress, 'update_interval') && ~isempty(progress.update_interval)
    interval = progress.update_interval;
end
end

function band = local_band(progress, field, index)
if isfield(progress, field) && numel(progress.(field)) >= index
    band = progress.(field)(index);
else
    band = index;
end
end
