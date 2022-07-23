% This function applies the boundary conditions to the linear system

function [K,R,fxnodes] = CalcBC(Tri,Kgel,Residual)

K=Kgel;
R=Residual;
fxnodes = ones(2*Tri(1).nnodes,1);          % d.o.f. (0:free ; 1:fixed/imposed)

%----AT BOTTOM--------AT BOTTOM--------AT BOTTOM--------AT BOTTOM----%
%----AT BOTTOM--------AT BOTTOM--------AT BOTTOM--------AT BOTTOM----%


%--------------------------------------------%
%   restraint x,y directions at the bottom   %
%                 (undercut)                 %
%--------------------------------------------%

% base = (2*(2*Tri(1).m+1));
% for i=1:base
% fxnodes(i,1) = 0;
%     for j=1:(2*Tri(1).nnodes)
%         K(i,j)=0;
%         if i==j
%             K(i,j)=1;
%             R(i,1)=0;
%         end
%     end
% end


%------------------------------------------------------%
%   restraint x (node 1) & y direction at the bottom   %
%                       (kerf)                         %
%------------------------------------------------------%

% base = ((2*Tri(1).m+1));
% for i=1:base
% fxnodes(i,1) = 0;
%     for j=1:(2*Tri(1).nnodes)
%         K(i,j) = mod(i,2)*K(i,j);
%         if i==j && mod(i,2)==0
%             K(i,j)=1;
%             R(i,1)=0;
%         end
%     end
% end
% 
% for j=1:(2*Tri(1).nnodes)
%     K(1,j) = 0;
% end
% K(1,1) = 1;
% R(1,1) = 0;

%----------------------------------------------------%
%   pontual restraint x,y directions at the bottom   %
%                      (dowel)                       %
%----------------------------------------------------%

% base = (2*(2*Tri(1).m+1));
% middle = (2*Tri(1).m+1);
% fxnodes(middle,1) = 0;
% fxnodes(middle+1,1) = 0;
% for j=1:(2*Tri(1).nnodes)
%     K(middle,j)=0;
%     K(middle+1,j)=0;
% end
% K(middle,middle)=1;
% K(middle+1,middle+1)=1;
% R(middle,1)=0;
% R(middle+1,1)=0;

%----AT THE TOP--------AT THE TOP--------AT THE TOP----%
%----AT THE TOP--------AT THE TOP--------AT THE TOP----%


%--------------------------------------%
%  restraint y directions on the top   %
%--------------------------------------%

% top = (2*(2*Tri(1).m+1));
% count = (2*Tri(1).n)*(2*(2*Tri(1).m+1));
% top = top + count;
% for i=(count+1):top
%     fxnodes(i,1) = mod(i,2);
%     for j=1:(2*Tri(1).nnodes)
%         K(i,j) = mod(i,2)*K(i,j);
%         if i==j && mod(i,2)==0
%             K(i,j)=1;
%             R(i,1)=0;
%         end
%     end
% end


%----------------------------------------%
%   isostatic slab fixed at the bottom   %
%----------------------------------------%

%base = 2*(2*Tri(1).m+1);
%middle = (2*Tri(1).m+1);
%middle = (2*Tri(1).m+1+2*2);
fx1 = 1;
fx2 = 2;
fx3 = 2*(2*Tri(1).m+1);
fxnodes(fx1,1) = 0;
fxnodes(fx2,1) = 0;
fxnodes(fx3,1) = 0;
for j=1:length(Tri(1).gnodes)
    K(fx1,j)=0;
    K(fx2,j)=0;
    K(fx3,j)=0;
end
K(fx1,fx1)=1;
K(fx2,fx2)=1;
K(fx3,fx3)=1;
R(fx1,1)=0;
R(fx2,1)=0;
R(fx3,1)=0;


end