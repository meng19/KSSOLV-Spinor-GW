function k_mu = page_solve_screened_w(coeff_pages, vmat)
%PAGE_SOLVE_SCREENED_W Solve reduced screened-W systems by frequency page.

nmu = size(coeff_pages, 1);
nfreq = size(coeff_pages, 3);
if nfreq == 1
    identity = eye(nmu, 'like', coeff_pages);
    k_mu = (identity - coeff_pages * vmat) \ coeff_pages;
    return;
end

if isa(coeff_pages, 'gpuArray')
    try
        system_pages = local_system_pages(coeff_pages, vmat);
        k_mu = pagefun(@mldivide, system_pages, coeff_pages);
        return;
    catch
    end
elseif exist('pagemldivide', 'file') == 2 || ...
        exist('pagemldivide', 'builtin') == 5
    try
        system_pages = local_system_pages(coeff_pages, vmat);
        k_mu = pagemldivide(system_pages, coeff_pages);
        return;
    catch
    end
end

k_mu = local_loop_solve(coeff_pages, vmat);
end

function system_pages = local_system_pages(coeff_pages, vmat)
nmu = size(coeff_pages, 1);
nfreq = size(coeff_pages, 3);
identity = eye(nmu, 'like', coeff_pages);

coeff_rows = reshape(permute(coeff_pages, [1, 3, 2]), ...
    nmu * nfreq, nmu);
coeff_v_rows = coeff_rows * vmat;
coeff_v_pages = permute(reshape(coeff_v_rows, nmu, nfreq, nmu), ...
    [1, 3, 2]);
system_pages = bsxfun(@plus, -coeff_v_pages, identity);
end

function k_mu = local_loop_solve(coeff_pages, vmat)
nmu = size(coeff_pages, 1);
nfreq = size(coeff_pages, 3);
k_mu = complex(zeros(nmu, nmu, nfreq, 'like', coeff_pages));
identity = eye(nmu, 'like', coeff_pages);
for ifreq = 1:nfreq
    coeff = coeff_pages(:, :, ifreq);
    system_matrix = identity - coeff * vmat;
    k_mu(:, :, ifreq) = system_matrix \ coeff;
end
end
