% It assembles the banded stiffness matrix

function K_bm = Assmbl_BM(delU,Tri,dN)
format long

%----------------%
%   BAND WIDTH   %
%----------------%

delta = 0;

for i=1:Tri(1).nel
    min = Tri(1).ConM(i,1);
    max = Tri(1).ConM(i,1);
    for j=1:6
        if min > Tri(1).ConM(i,j)
            min = Tri(1).ConM(i,j);
        end
        if max < Tri(1).ConM(i,j)
            max = Tri(1).ConM(i,j);
        end
    end
    if delta < (max-min)
        delta = max - min;
    end
end
band_w = delta+1;                         % half band width = ma + 1
ma = 2*band_w - 1;



%-----------------------------------------%
%   ASSEMBLY OF BANDED STIFFNESS MATRIX   %
%-----------------------------------------%

K_bm = zeros(2*Tri(1).nnodes,(2*ma+1));        % declaration of banded stiffness matrix in global coordinates

for count=1:Tri(1).nel
    %K = CalcStiff_WI(Tri,dN,count);                % element stiffness matrix in local coordinates
    K = CalcStiff_NL(delU,Tri,dN,count);
    for node_i=1:6
        iglobal = Tri(1).ConM(count,node_i);    % i_bm = iglobal
        for node_j=1:6
            jglobal = Tri(1).ConM(count,node_j);
            j_bm = ma +1 + (2*jglobal-1) - (2*iglobal-1);
            K_bm((2*iglobal-1),(j_bm)) = K((2*node_i-1),(2*node_j-1)) + K_bm((2*iglobal-1),(j_bm));
            K_bm((2*iglobal-1),(j_bm+1)) = K((2*node_i-1),(2*node_j)) + K_bm((2*iglobal-1),(j_bm+1));
            j_bm = j_bm - 1;
            K_bm(2*iglobal,(j_bm)) = K(2*node_i,(2*node_j-1)) + K_bm(2*iglobal,(j_bm));
            K_bm(2*iglobal,j_bm+1) = K(2*node_i,2*node_j) + K_bm(2*iglobal,j_bm+1);
        end
    end
end
end