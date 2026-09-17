function ind_mu = partial_qrcp_mex_sample(products, rank_mu)
%PARTIAL_QRCP_MEX_SAMPLE Select points with the compiled truncated QRCP.

if exist('partial_qrcp_mex', 'file') ~= 3
    error('ISDF:PartialQRCPMexMissing', ...
        ['partial_qrcp_mex is not built. Run ' ...
         'isdf.build_partial_qrcp_mex before using this sample method.']);
end
ind_mu = partial_qrcp_mex(gather_if_gpu(products), rank_mu);
end
