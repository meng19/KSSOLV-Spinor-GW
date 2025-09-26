function [sys, options, syms] = read_kssolv_gw(sys, H0, X0)
%% 
% For now, there is no symmetry in kssolv
syms.ntran = 1;
syms.ntranq = 0;
syms.mtrx = {[1,0,0;0,1,0;0,0,1];[1,0,0;0,-1,0;0,0,-1];[1,0,0;0,0,1;0,-1,0];[1,0,0;0,0,-1;0,1,0];[1,0,0;0,1,0;0,0,-1];[1,0,0;0,-1,0;0,0,1];[1,0,0;0,0,-1;0,-1,0];[1,0,0;0,0,1;0,1,0];[0,0,1;0,-1,0;1,0,0];[0,0,-1;0,-1,0;-1,0,0];[0,0,1;0,1,0;-1,0,0];[0,0,-1;0,1,0;1,0,0];[-1,0,0;0,0,1;0,1,0];[-1,0,0;0,0,-1;0,-1,0];[-1,0,0;0,1,0;0,0,-1];[-1,0,0;0,-1,0;0,0,1];[0,1,0;0,0,1;1,0,0];[0,-1,0;0,0,1;-1,0,0];[0,1,0;0,0,-1;-1,0,0];[0,-1,0;0,0,-1;1,0,0];[0,0,1;1,0,0;0,1,0];[0,0,1;-1,0,0;0,-1,0];[0,0,-1;-1,0,0;0,1,0];[0,0,-1;1,0,0;0,-1,0];[-1,0,0;0,-1,0;0,0,-1];[0,1,0;1,0,0;0,0,-1];[0,-1,0;-1,0,0;0,0,-1];[-1,0,0;0,1,0;0,0,1];[0,-1,0;-1,0,0;0,0,1];[0,1,0;1,0,0;0,0,1];[0,-1,0;1,0,0;0,0,-1];[0,1,0;-1,0,0;0,0,-1];[0,0,-1;0,1,0;-1,0,0];[0,0,1;0,1,0;1,0,0];[0,0,-1;0,-1,0;1,0,0];[0,0,1;0,-1,0;-1,0,0];[0,1,0;-1,0,0;0,0,1];[0,-1,0;1,0,0;0,0,1];[-1,0,0;0,0,-1;0,1,0];[-1,0,0;0,0,1;0,-1,0];[0,-1,0;0,0,-1;-1,0,0];[0,1,0;0,0,-1;1,0,0];[0,-1,0;0,0,1;1,0,0];[0,1,0;0,0,1;-1,0,0];[0,0,-1;-1,0,0;0,-1,0];[0,0,-1;1,0,0;0,1,0];[0,0,1;1,0,0;0,-1,0];[0,0,1;-1,0,0;0,1,0]};
syms.nrot = 1;
syms.indsub = zeros(48,1);
syms.kgzero = zeros(48,3);

%% ev
% Unit Ry
num_vectors = size(H0.eband, 1);
vector_length = length(H0.eband{1,1});
matrix = zeros(vector_length, num_vectors);

for i = 1:num_vectors
    matrix(:, i) = H0.eband{i,1};
end
options.ev = 2 * matrix;

%% 
options.efermi = 2 * sys.efermi;
options.rho0 = H0.rho;
options.X0 = X0;

%% 
for ik = 1:sys.nkpts
    for ispin = 1:sys.nspin
        occ(:, ispin) = X0.wavefuncell{ik, ispin}.occ;
        ifmax = find(occ(:, ispin) >= 0.5, 1, 'last' );
        options.ifmax(ik, ispin) = ifmax;
    end
end

% For now, there is only one k=0 G-grid in kssolv
sigrid = Ggrid(sys, sys.ecut2);
mill = sigrid.kkxyz(X0.wavefuncell{1, 1}.idxnz, :);
mill = round(mill / sys.bvec, 6); % reduce accuracy to prevent truncation errors
for ik = 1:sys.nkpts
    options.mill{1, ik} = mill;
end

%% Calculate vxc of every band from rho and wav
vxc = cal_vxc(sys, options);
sys.vxc = vxc;

%% Others
if (sys.lsda == 1)
    nspin = 2;
    if (sys.lspinorb == 1)
        sys.nspinor = 2;
    else
        sys.nspinor = 1;
    end
else
    nspin = 1;
    if (sys.lspinorb == 1)
        sys.nspinor = 2;
    else
        sys.nspinor = 1;
    end
end

options.kpts = round(sys.kpts / sys.bvec, 6); % reduce accuracy to prevent truncation errors
end
