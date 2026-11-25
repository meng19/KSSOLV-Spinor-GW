function sig = sigma_set_defaults(sig)
% ==================== 设置 sig 的默认值 ====================
if ~isfield(sig, 'use_gpu')
    sig.use_gpu = false;                % 默认不使用GPU
end

if ~isfield(sig, 'precompute_wav')
    sig.precompute_wav = false;
end

% Optional epsilon checkpoint used when sigma is run in a later session.
if ~isfield(sig, 'epsilon_checkpoint_file') || isempty(sig.epsilon_checkpoint_file)
    sig.epsilon_checkpoint_file = '';
end
if isfield(sig, 'epsilon_checkpoint') && ~isempty(sig.epsilon_checkpoint)
    sig.epsilon_checkpoint_file = sig.epsilon_checkpoint;
end
if ~(ischar(sig.epsilon_checkpoint_file) || ...
        (isstring(sig.epsilon_checkpoint_file) && ...
        isscalar(sig.epsilon_checkpoint_file)))
    error('GW:EpsilonCheckpointFile', ...
        'sig.epsilon_checkpoint_file must be a character vector or string scalar.');
end

if ~isfield(sig, 'freq_dep')
    sig.freq_dep = 0;                  % 默认静态COHSEX
end

if ~isfield(sig, 'freq_dep_method')
    sig.freq_dep_method = 2;           % 默认方法
end

if ~isfield(sig, 'no_symmetries_q_grid')
    sig.no_symmetries_q_grid = false;
end

if ~isfield(sig, 'exact_static_ch')
    sig.exact_static_ch = false;
end

if ~isfield(sig, 'isdf') || isempty(sig.isdf) || ~isstruct(sig.isdf)
    sig.isdf = struct();
end
sig.isdf = isdf.workflow_defaults(sig.isdf, sig.freq_dep);

% The legacy path uses the NN product space for every sigma contribution.
% A VN space can be selected explicitly for the bare-exchange term only.
if ~isfield(sig.isdf, 'exchange_space') || isempty(sig.isdf.exchange_space)
    sig.isdf.exchange_space = 'nn';
end
if ~isfield(sig.isdf, 'global_nn_space') || isempty(sig.isdf.global_nn_space)
    sig.isdf.global_nn_space = false;
end
if ~isfield(sig.isdf, 'global_vn_space') || isempty(sig.isdf.global_vn_space)
    sig.isdf.global_vn_space = false;
end
if ~isfield(sig.isdf, 'reuse_nn_for_vn') || isempty(sig.isdf.reuse_nn_for_vn)
    sig.isdf.reuse_nn_for_vn = false;
end
if ~isfield(sig.isdf, 'reuse_eps_real_wfn') || isempty(sig.isdf.reuse_eps_real_wfn)
    sig.isdf.reuse_eps_real_wfn = false;
end
if ~isfield(sig.isdf, 'right_band_block_size') || ...
        isempty(sig.isdf.right_band_block_size)
    sig.isdf.right_band_block_size = Inf;
elseif ischar(sig.isdf.right_band_block_size) || ...
        (isstring(sig.isdf.right_band_block_size) && ...
        isscalar(sig.isdf.right_band_block_size))
    sig.isdf.right_band_block_size = str2double( ...
        char(sig.isdf.right_band_block_size));
end
if ~(isnumeric(sig.isdf.right_band_block_size) && ...
        isscalar(sig.isdf.right_band_block_size) && ...
        sig.isdf.right_band_block_size >= 1)
    error('ISDF:RightBandBlockSize', ...
        'sig.isdf.right_band_block_size must be a positive scalar or Inf.');
end
if ~isfield(sig.isdf, 'component_cache_gb') || ...
        isempty(sig.isdf.component_cache_gb)
    sig.isdf.component_cache_gb = Inf;
end
if ~(isnumeric(sig.isdf.component_cache_gb) && ...
        isscalar(sig.isdf.component_cache_gb) && ...
        sig.isdf.component_cache_gb >= 0)
    error('ISDF:ComponentCacheLimit', ...
        'sig.isdf.component_cache_gb must be a non-negative scalar or Inf.');
end
if ~isfield(sig.isdf, 'fullfreq_batch_bytes') || ...
        isempty(sig.isdf.fullfreq_batch_bytes)
    % Memory budget for one fullfreq contraction batch, which holds
    % rank^2 x batch_size complex coefficients.
    sig.isdf.fullfreq_batch_bytes = 48 * (1024^3);
else
    sig.isdf.fullfreq_batch_bytes = local_parse_byte_size( ...
        sig.isdf.fullfreq_batch_bytes, 'sig.isdf.fullfreq_batch_bytes');
end
if ~(isnumeric(sig.isdf.fullfreq_batch_bytes) && ...
        isscalar(sig.isdf.fullfreq_batch_bytes) && ...
        sig.isdf.fullfreq_batch_bytes >= 16 * 16)
    error('ISDF:FullfreqBatchBytes', ...
        ['sig.isdf.fullfreq_batch_bytes must be a positive byte size, ' ...
         'such as "48GB".']);
end

if ~isfield(sig.isdf, 'screened_kernel_cache_gb') || ...
        isempty(sig.isdf.screened_kernel_cache_gb)
    % Cache budget for reusing reduced screened kernels shared by several
    % diagonal bands.  Inf keeps every reusable kernel (no memory cap); a
    % non-positive value disables the reuse.  The value is in GB.
    if isfield(sig.isdf, 'screened_kernel_cache_bytes') && ...
            ~isempty(sig.isdf.screened_kernel_cache_bytes)
        sig.isdf.screened_kernel_cache_gb = ...
            sig.isdf.screened_kernel_cache_bytes / (1024^3);
    else
        sig.isdf.screened_kernel_cache_gb = Inf;
    end
end
if isfield(sig.isdf, 'screened_kernel_cache_bytes')
    % Backward compatibility: keep legacy alias in bytes, but prefer the GB
    % field for user input.
    sig.isdf.screened_kernel_cache_bytes = ...
        sig.isdf.screened_kernel_cache_gb * (1024^3);
end
if ~(isnumeric(sig.isdf.screened_kernel_cache_gb) && ...
        isscalar(sig.isdf.screened_kernel_cache_gb) && ...
        sig.isdf.screened_kernel_cache_gb >= 0)
    error('ISDF:ScreenedKernelCacheLimit', ...
        ['sig.isdf.screened_kernel_cache_gb must be a ' ...
         'non-negative scalar or Inf.']);
end
end

function value = local_parse_byte_size(value, option_name)
% Parse byte sizes with decimal GB or binary GiB style units.
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    return;
end
text = lower(strtrim(char(value)));
if strcmp(text, 'inf')
    value = Inf;
    return;
end
numeric_value = str2double(text);
if ~isnan(numeric_value)
    value = numeric_value;
else
    number_pattern = ...
        '([+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:e[+-]?[0-9]+)?)';
    tokens = regexp(text, ...
        ['^' number_pattern '\s*([kmgt]?)(i?)b$'], 'tokens', 'once');
    if isempty(tokens) || (isempty(tokens{2}) && ~isempty(tokens{3}))
        error('ISDF:SizeOption', ...
            ['%s must be a positive byte size, such as "48GB", ' ...
            '"48GiB", or "4.8e10".'], option_name);
    end
    value = str2double(tokens{1});
    unit = tokens{2};
    if isempty(unit)
        multiplier = 1;
    else
        if isempty(tokens{3})
            base = 1000;
        else
            base = 1024;
        end
        switch unit
            case 'k'
                exponent = 1;
            case 'm'
                exponent = 2;
            case 'g'
                exponent = 3;
            case 't'
                exponent = 4;
        end
        multiplier = base^exponent;
    end
    value = value * multiplier;
end
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        ~isnan(value) && value >= 1)
    error('ISDF:SizeOption', '%s must be a positive byte size.', option_name);
end
end
