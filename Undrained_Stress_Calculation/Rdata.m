% This fucntion read input data (number of nodes, structure length & width)

function [Tri] = Rdata(m,n)

Tri(1).b = 0.03;                                           %(m)
Tri(1).h = 1.03;                                           %(m)
%m = 5;                                             % row subdivision
%n = 50;                                            % column subdivision
Tri(1).nel = 2*m*n;                                 % number of elements
Tri(1).nnodes = (9+(m-1)*6)*n-(2*m+1)*(n-1);        % number of nodes
Tri(1).ConM = zeros(Tri(1).nel,6);                  % connectivity matrix (it connects local to global nodes)
Tri(1).gnodes_0 = zeros(2*Tri(1).nnodes,1);         % global nodes coordinates in reference system
Tri(1).gnodes = zeros(2*Tri(1).nnodes,1);           % global nodes coordinates in updated system

%dcalcite = 125*10^-6;       %average size of calcite grain in meter (Ferrero et al 2014)

delm = Tri(1).b/m;
deln = Tri(1).h/n;
nlin = (2*m+1); % number of nodes in each line
ncol = (2*n+1); % number of nodes in each column

%--------------------------------%
%   NODAL COORDINATES (GLOBAL)   %
%--------------------------------%

temp=1;
for j=1:n
    for k=1:2
        if k==1
            Tri(1).gnodes_0(temp,1) = 0;                 % direction "X"
            Tri(1).gnodes_0(temp+1,1) = (j-1)*deln;      % direction "Y"
            temp=temp+2;
            for i=1:m
                Tri(1).gnodes_0(temp,1) = (i-0.5)*delm;  % direction "X"
                Tri(1).gnodes_0(temp+1,1) = (j-1)*deln;  % direction "Y"
                temp=temp+2;
                Tri(1).gnodes_0(temp,1) = (i)*delm;      % direction "X"
                Tri(1).gnodes_0(temp+1,1) = (j-1)*deln;  % direction "Y"
                temp=temp+2 ;
            end
        else
            Tri(1).gnodes_0(temp,1) = 0;                 % direction "X"
            Tri(1).gnodes_0(temp+1,1) = (j-0.5)*deln;    % direction "Y"
            temp=temp+2;
            for i=1:(m)
                Tri(1).gnodes_0(temp,1) = (i-0.5)*delm;  % direction "X"
                Tri(1).gnodes_0(temp+1,1) = (j-0.5)*deln;% direction "Y"
                temp=temp+2;
                Tri(1).gnodes_0(temp,1) = (i)*delm;      % direction "X"
                Tri(1).gnodes_0(temp+1,1) = (j-0.5)*deln;% direction "Y"
                temp=temp+2;
            end
        end
    end
end
j=n;
Tri(1).gnodes_0(temp,1) = 0;                 % direction "X"
Tri(1).gnodes_0(temp+1,1) = (j)*deln;        % direction "Y"
temp=temp+2;
for i=1:(m)
    Tri(1).gnodes_0(temp,1) = (i-0.5)*delm;  % direction "X"
    Tri(1).gnodes_0(temp+1,1) = (j)*deln;    % direction "Y"
    temp=temp+2;
    Tri(1).gnodes_0(temp,1) = (i)*delm;      % direction "X"
    Tri(1).gnodes_0(temp+1,1) = (j)*deln;    % direction "Y"
    temp=temp+2;
end

%-------------------------------%
%   NODAL COORDINATES (LOCAL)   %
%-------------------------------%

count=0;
for j=1:n
    if mod(j,2)==1
        for i=1:m
            if mod(i,2)==1
                count = count + 1;
                [Tri,count] = elem34(Tri,i,j,m,n,count);
            else
                count = count + 1;
                [Tri,count] = elem12(Tri,i,j,m,n,count);
            end
        end
    else
        for i=1:m
            if mod(i,2)==0
                count = count + 1;
                [Tri,count] = elem34(Tri,i,j,m,n,count);
            else
                count = count + 1;
                [Tri,count] = elem12(Tri,i,j,m,n,count);
            end
        end
    end
end

Tri(1).gnodes = Tri(1).gnodes_0;


%--------------------------%
%   PLOT TRIANGULAR MESH   %
%--------------------------%

for i=1:Tri(1).nel
    for j=1:3
        temp = Tri(1).ConM(i,j);
        T1(i,j) = temp;
        xx(temp) = Tri(1).gnodes_0(2*temp-1);
        yy(temp) = Tri(1).gnodes_0(2*temp);
    end
end
triplot(T1,xx,yy)
axis([0, 0.4, 0, 0.4]);
axis([0, 0.4, 0, 0.4])
end