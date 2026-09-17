function build_partial_qrcp_mex
%BUILD_PARTIAL_QRCP_MEX Build the blocked, truncated QRCP MEX kernel.

source_dir = fullfile(fileparts(mfilename('fullpath')), 'private');
build_mex_kernel(fullfile(source_dir, 'partial_qrcp_mex.cpp'), source_dir);
end
