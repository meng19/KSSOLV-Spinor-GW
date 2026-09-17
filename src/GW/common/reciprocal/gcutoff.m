function [gcutoff] = gcutoff(ng, ekin, isrtrq, ecutoff, method)
%GCUTOFF Number of G vectors whose kinetic energy stays within the cutoff.
%   EKIN(ISRTRQ) is sorted ascending (SORTRX orders by kinetic energy
%   first), so the cutoff index is the last position whose sorted energy
%   does not exceed ECUTOFF.  ISRTRQ stays in the signature for caller
%   compatibility.
%
%   Two interchangeable implementations are kept; switch between them per
%   call with GCUTOFF(..., METHOD), or globally by editing DEFAULT_METHOD
%   below:
%     'binary' (default)  Upper-midpoint binary search, O(log(ng)) probes
%                         into EKIN(ISRTRQ).  Flat microsecond cost; best
%                         for the large full-grid sizes used in GW.  (The
%                         legacy round((gup+gdn)/2) midpoint equals GUP
%                         once the interval held two entries, which made
%                         the loop spin for O(ng) iterations.)
%     'count'             Single vectorized pass NNZ(EKIN <= ECUTOFF).
%                         Best for small grids or EKIN already on the GPU.

DEFAULT_METHOD = 'binary';
if nargin >= 5 && ~isempty(method)
    selected = lower(char(method));
else
    selected = DEFAULT_METHOD;
end

switch selected
    case 'binary'
        lo = 1;
        hi = ng;
        while lo < hi
            mid = floor((lo + hi + 1) / 2);
            if ekin(isrtrq(mid)) <= ecutoff
                lo = mid;
            else
                hi = mid - 1;
            end
        end
        gcutoff = lo;
    case 'count'
        gcutoff = nnz(ekin(:) <= ecutoff);
        if gcutoff < 1
            % Legacy behavior: even when no vector passes the cutoff the
            % first (smallest-energy) G vector is always retained.
            gcutoff = 1;
        end
    otherwise
        error('GW:UnknownGcutoffMethod', ...
            'gcutoff method must be ''binary'' or ''count'' (got ''%s'').', ...
            selected);
end
end
