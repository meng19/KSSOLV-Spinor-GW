function build_mex_kernel(source_file, output_dir, varargin)
%BUILD_MEX_KERNEL Build one C/C++ MEX kernel with consistent options.
%   BUILD_MEX_KERNEL(SOURCE_FILE, OUTPUT_DIR) clears a stale same-platform
%   binary and invokes MEX using the interleaved-complex API and optimization.
%   Extra arguments are passed through to MEX, for example '-lmwblas'.

if nargin < 2 || isempty(output_dir)
    output_dir = fileparts(source_file);
end
if ~isfile(source_file)
    error('KSSOLV:MexSourceMissing', 'MEX source does not exist: %s', source_file);
end

[~, stem] = fileparts(source_file);
stale = fullfile(output_dir, [stem '.', mexext]);
if isfile(stale)
    clear(stem);
    delete(stale);
end
mex('-R2018a', '-O', '-outdir', output_dir, source_file, varargin{:});
end
