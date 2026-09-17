// PARTIAL_QRCP_MEX  Rank-truncated blocked QRCP using LAPACK ZLAQPS.
// Build with build_partial_qrcp_mex.  The input is products (ngrid-by-nsketch)
// and is interpreted as products.' without materializing that transpose.
// Factorization stops after rank_mu pivots, using ZGEQP3's LAPACK panel size.

#include "mex.h"
#include "blas.h"
#include "lapack.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

namespace {

using Complex = std::complex<double>;

double squared_norm(const std::vector<Complex>& a, mwSize m, mwSize column,
                    mwSize first_row) {
    double value = 0.0;
    const mwSize offset = column * m;
    for (mwSize row = first_row; row < m; ++row) {
        value += std::norm(a[offset + row]);
    }
    return value;
}

void swap_columns(std::vector<Complex>& a, mwSize m, mwSize first,
                  mwSize second) {
    if (first == second) {
        return;
    }
    const mwSize first_offset = first * m;
    const mwSize second_offset = second * m;
    for (mwSize row = 0; row < m; ++row) {
        std::swap(a[first_offset + row], a[second_offset + row]);
    }
}

}  // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs,
                 const mxArray* prhs[]) {
    if (nrhs != 2 || nlhs > 1) {
        mexErrMsgIdAndTxt("ISDF:PartialQRCPMexUsage",
            "Usage: pivots = partial_qrcp_mex(products, rank_mu).");
    }
    const mxArray* products = prhs[0];
    if (!mxIsDouble(products) || mxIsSparse(products) || mxGetNumberOfDimensions(products) != 2) {
        mexErrMsgIdAndTxt("ISDF:PartialQRCPMexInput",
            "products must be a full double 2-D matrix.");
    }
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
            mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("ISDF:PartialQRCPMexRank",
            "rank_mu must be a real scalar.");
    }

    const mwSize ngrid = mxGetM(products);
    const mwSize nsketch = mxGetN(products);
    const double requested_rank = mxGetScalar(prhs[1]);
    if (!std::isfinite(requested_rank) || requested_rank < 1 ||
            std::floor(requested_rank) != requested_rank) {
        mexErrMsgIdAndTxt("ISDF:PartialQRCPMexRank",
            "rank_mu must be a positive integer.");
    }
    const mwSize rank_mu = static_cast<mwSize>(requested_rank);
    if (rank_mu > std::min(ngrid, nsketch)) {
        mexErrMsgIdAndTxt("ISDF:PartialQRCPMexRank",
            "rank_mu cannot exceed min(size(products)).");
    }

    // MATLAB stores products(j, i) at j + ngrid*i, exactly the storage
    // position of products.'(i, j) when interpreted as nsketch-by-ngrid.
    const mwSize count = ngrid * nsketch;
    std::vector<Complex> a(count);
    if (mxIsComplex(products)) {
        const mxComplexDouble* input = mxGetComplexDoubles(products);
        for (mwSize index = 0; index < count; ++index) {
            a[index] = Complex(input[index].real, input[index].imag);
        }
    } else {
        const double* input = mxGetDoubles(products);
        for (mwSize index = 0; index < count; ++index) {
            a[index] = Complex(input[index], 0.0);
        }
    }

    const ptrdiff_t m = static_cast<ptrdiff_t>(nsketch);
    const ptrdiff_t n = static_cast<ptrdiff_t>(ngrid);
    const ptrdiff_t increment = 1;
    std::vector<ptrdiff_t> pivots(ngrid);
    std::vector<double> vn1(ngrid);
    std::vector<double> vn2(ngrid);
    std::vector<Complex> tau(rank_mu);
    for (mwSize column = 0; column < ngrid; ++column) {
        pivots[column] = static_cast<ptrdiff_t>(column + 1);
        vn1[column] = dznrm2(&m, reinterpret_cast<const double*>(
            a.data() + column * nsketch), &increment);
        vn2[column] = vn1[column];
    }

    // Match ZGEQP3: NB = ILAENV(1, 'ZGEQRF', ' ', M, N, -1, -1).
    const ptrdiff_t query = 1;
    const ptrdiff_t minus_one = -1;
    const char routine[] = "ZGEQRF";
    const char options[] = " ";
    const ptrdiff_t block_size = std::max<ptrdiff_t>(1, ilaenv(
        &query, routine, options, &m, &n, &minus_one, &minus_one, 6, 1));
    ptrdiff_t start = 0;
    while (start < static_cast<ptrdiff_t>(rank_mu)) {
        const ptrdiff_t remaining_columns = n - start;
        const ptrdiff_t panel_width = std::min(block_size,
            static_cast<ptrdiff_t>(rank_mu) - start);
        ptrdiff_t factored = 0;
        std::vector<Complex> auxiliary(panel_width);
        std::vector<Complex> f(remaining_columns * panel_width);
        zlaqps(&m, &remaining_columns, &start, &panel_width, &factored,
            reinterpret_cast<double*>(a.data() + start * nsketch), &m,
            pivots.data() + start, reinterpret_cast<double*>(tau.data() + start),
            vn1.data() + start, vn2.data() + start,
            reinterpret_cast<double*>(auxiliary.data()),
            reinterpret_cast<double*>(f.data()), &remaining_columns);
        if (factored <= 0) {
            mexErrMsgIdAndTxt("ISDF:PartialQRCPMexRankDeficient",
                "The QRCP sketch has rank below the requested ISDF rank.");
        }
        start += factored;
    }

    plhs[0] = mxCreateDoubleMatrix(1, rank_mu, mxREAL);
    double* output = mxGetDoubles(plhs[0]);
    for (mwSize index = 0; index < rank_mu; ++index) {
        output[index] = static_cast<double>(pivots[index]);
    }
}
