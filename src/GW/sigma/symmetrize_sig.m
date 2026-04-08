function sig = symmetrize_sig(sig, emf, ndiag_min, ndiag_max, nspin, fieldNames, tolerance)
% 默认参数
if nargin < 7, tolerance = 1e-5; end
if nargin < 6 || isempty(fieldNames), fieldNames = {'cor', 'sig'}; end

for ik = 1:sig.nkn
    for ispin = 1:nspin
        % 找出简并分组
        nsubs = 1;
        ndeg = [];
        ndeg(nsubs) = 1;
        
        for ii = ndiag_min+1 : ndiag_max
            if abs(emf(ii,ik,ispin) - emf(ii-1,ik,ispin)) < tolerance
                ndeg(nsubs) = ndeg(nsubs) + 1;
            else
                nsubs = nsubs + 1;
                ndeg(nsubs) = 1;
            end
        end
        
        % 对每个简并组求平均
        istop = 0;
        for ii = 1:nsubs
            istart = istop + 1;
            istop = istart + ndeg(ii) - 1;
            
            % 对所有需要对称化的场求平均
            for iField = 1:length(fieldNames)
                fname = fieldNames{iField};
                avg = mean(sig.(fname)(istart:istop, ik, ispin));
                sig.(fname)(istart:istop, ik, ispin) = avg;
            end
        end
    end
end
end