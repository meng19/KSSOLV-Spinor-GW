function eps = epsilon_warn_cauchy_fallback(ctx, eps)
%EPSILON_WARN_CAUCHY_FALLBACK Summarize reduced Cauchy fallback pages.

if ~strcmp(ctx.method, 'reduced_basis')
    return;
end
if ~isfield(ctx.eps, 'isdf') || ...
        ~strcmpi(ctx.eps.isdf.reduced_solver, 'cauchy')
    return;
end
if ~isfield(eps, 'isdf_reduced_info') || isempty(eps.isdf_reduced_info)
    return;
end

summary = repmat(local_empty_summary(), ctx.nq, 1);
warning_lines = {};
for iq = 1:ctx.nq
    nfreq = ctx.pol.nfreq;
    nreal = local_real_frequency_count(ctx.pol, nfreq);
    qdata = ctx.qdata{iq};
    fallback_pages = false(1, nfreq);
    block_counts = zeros(ctx.nspin, qdata.nrq);
    block_real_counts = zeros(ctx.nspin, qdata.nrq);

    for ispin = 1:ctx.nspin
        for ik = 1:qdata.nrq
            info = eps.isdf_reduced_info{iq, ispin, ik};
            pages = local_fallback_pages(info, nfreq);
            if isempty(pages)
                continue;
            end
            fallback_pages = fallback_pages | pages;
            block_counts(ispin, ik) = nnz(pages);
            block_real_counts(ispin, ik) = nnz(pages(1:nreal));
        end
    end

    real_pages = fallback_pages(1:nreal);
    summary(iq).k_index = iq;
    summary(iq).q_vector = ctx.pol.qpt(iq, :);
    summary(iq).real_pages = real_pages;
    summary(iq).real_count = nnz(real_pages);
    summary(iq).all_pages = fallback_pages;
    summary(iq).all_count = nnz(fallback_pages);
    summary(iq).block_real_counts = block_real_counts;
    summary(iq).block_counts = block_counts;
    summary(iq).affected_blocks = nnz(block_counts);
    summary(iq).total_blocks = numel(block_counts);

    if summary(iq).all_count > 0
        warning_lines{end + 1} = local_warning_line( ...
            summary(iq), nreal, nfreq); %#ok<AGROW>
    end
end

eps.isdf_cauchy_fallback_summary = summary;
if ~isempty(warning_lines)
    local_warn_summary(warning_lines);
end
end

function summary = local_empty_summary()
summary = struct( ...
    'k_index', [], ...
    'q_vector', [], ...
    'real_pages', [], ...
    'real_count', 0, ...
    'all_pages', [], ...
    'all_count', 0, ...
    'block_real_counts', [], ...
    'block_counts', [], ...
    'affected_blocks', 0, ...
    'total_blocks', 0);
end

function nreal = local_real_frequency_count(pol, nfreq)
if isfield(pol, 'nfreq_rel') && ~isempty(pol.nfreq_rel)
    nreal = min(pol.nfreq_rel, nfreq);
else
    nreal = nfreq;
end
end

function pages = local_fallback_pages(info, nfreq)
pages = [];
if ~isstruct(info) || ~isfield(info, 'fallback_pages') || ...
        isempty(info.fallback_pages)
    return;
end
pages = logical(info.fallback_pages(:).');
if numel(pages) < nfreq
    pages(end + 1:nfreq) = false;
elseif numel(pages) > nfreq
    pages = pages(1:nfreq);
end
end

function line = local_warning_line(summary, nreal, nfreq)
qvec = summary.q_vector;
line = sprintf( ...
    'K %d q=(%.4f, %.4f, %.4f): real=%d/%d page(s), blocks=%d/%d', ...
    summary.k_index, qvec(1), qvec(2), qvec(3), ...
    summary.real_count, nreal, ...
    summary.affected_blocks, summary.total_blocks);
if summary.all_count ~= summary.real_count || nfreq ~= nreal
    line = sprintf('%s, all=%d/%d page(s)', ...
        line, summary.all_count, nfreq);
end
end

function local_warn_summary(warning_lines)
backtrace_state = warning('query', 'backtrace');
cleanup = onCleanup( ...
    @() warning(backtrace_state.state, 'backtrace')); %#ok<NASGU>
warning('off', 'backtrace');
message = sprintf('%s\n', warning_lines{:});
warning('ISDF:CauchyFrequencyDirectFallback', ...
    ['Cauchy contour fallback summary by K/q block:' newline ...
     '%s' ...
     'Direct reduced summation was used only for the listed page(s).'], ...
    message);
end
