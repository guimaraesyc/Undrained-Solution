% This function applies the boundary conditions to the linear system

function [K,R,fxnodes] = CalcBC_BM(Tri,K_bm,Residual)

K = K_bm;
R = Residual;
fxnodes = ones(2*Tri(1).nnodes,1);          % d.o.f. (1:free ; 0:fixed/imposed)
[row,col] = size(K);
ma = 0.5*(col+1) - 1;

%----AT BOTTOM--------AT BOTTOM--------AT BOTTOM--------AT BOTTOM----%
%----AT BOTTOM--------AT BOTTOM--------AT BOTTOM--------AT BOTTOM----%


%--------------------------------------------%
%   restraint x,y directions at the bottom   %
%                 (undercut)                 %
%--------------------------------------------%

% base = (2*(2*Tri(1).m+1));
% for i=1:base
% fxnodes(i,1) = 0;
%     for j=1:(2*ma+1)
%         K(i,j)=0;
%         if j==(ma+1)
%             K(i,ma+1)=1;
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
%     for j=1:(2*ma+1)
%         K(i,j) = mod(i,2)*K(i,j);
%         if j==(ma+1) && mod(i,2)==0
%             K(i,j)=1;
%             R(i,1)=0;
%         end
%     end
% end
% 
% for j=1:(2*ma+1)
%     K(1,j) = 0;
% end
% K(1,ma+1) = 1;
% R(1,1) = 0;

%----------------------------------------------------%
%   pontual restraint x,y directions at the bottom   %
%               (EXCENTRIC dowel)                    %
%----------------------------------------------------%

% base = (2*(2*Tri(1).m+1));
% middle = (2*Tri(1).m+1+8);
% %middle = (2*Tri(1).m+1+2*2);
% fxnodes(middle,1) = 0;
% fxnodes(middle+1,1) = 0;
% for j=1:(2*ma+1)
%     K(middle,j)=0;
%     K(middle+1,j)=0;
% end
% K(middle,ma+1)=1;
% K(middle+1,ma+1)=1;
% R(middle,1)=0;
% R(middle+1,1)=0;
% 
% base = (2*(2*Tri(1).m+1));
% middle = (2*Tri(1).m+1);
% fxnodes(middle,1) = 0;
% fxnodes(middle+1,1) = 0;
% for j=1:(2*ma+1)
%     K(middle,j)=0;
%     K(middle+1,j)=0;
% end
% K(middle,ma+1)=1;
% K(middle+1,ma+1)=1;
% R(middle,1)=0;
% R(middle+1,1)=0;

%----------------------------------------------------%
%   pontual restraint x,y directions at the bottom   %
%                      (dowel)                       %
%----------------------------------------------------%

% base = (2*(2*Tri(1).m+1));
% middle = (2*Tri(1).m+1);
% %middle = (2*Tri(1).m+1+2*2);
% fxnodes(middle,1) = 0;
% fxnodes(middle+1,1) = 0;
% for j=1:(2*ma+1)
%     K(middle,j)=0;
%     K(middle+1,j)=0;
% end
% K(middle,ma+1)=1;
% K(middle+1,ma+1)=1;
% R(middle,1)=0;
% R(middle+1,1)=0;
% 
% base = (2*(2*Tri(1).m+1));
% middle = (2*Tri(1).m+1);
% fxnodes(middle,1) = 0;
% fxnodes(middle+1,1) = 0;
% for j=1:(2*ma+1)
%     K(middle,j)=0;
%     K(middle+1,j)=0;
% end
% K(middle,ma+1)=1;
% K(middle+1,ma+1)=1;
% R(middle,1)=0;
% R(middle+1,1)=0;

% %----AT THE TOP--------AT THE TOP--------AT THE TOP----%
% %----AT THE TOP--------AT THE TOP--------AT THE TOP----%
% 
% 
% %--------------------------------------%
% %  restraint y directions at the top   %
% %--------------------------------------%

% top = (2*(2*Tri(1).m+1));
% temp = (2*Tri(1).n)*top;
% top = top + temp;
% for i=(temp+1):top
%     fxnodes(i,1) = mod(i,2);
%     for j=1:(2*ma+1)
%         K(i,j) = mod(i,2)*K(i,j);
%         if j==(ma+1) && mod(i,2)==0
%             K(i,ma+1)=1;
%             R(i,1)=0;
%         end
%     end
% end



%----------------------------------------%
%   isostatic slab fixed at the bottom   %
%----------------------------------------%

% fx1 = 1;
% fx2 = 2;
% fx3 = (2*(2*Tri(1).m+1));
% fxnodes(fx1,1) = 0;
% fxnodes(fx2,1) = 0;
% fxnodes(fx3,1) = 0;
% 
% for j=1:(2*ma+1)
%     K(fx1,j) = 0;
%     K(fx2,j) = 0;
%     K(fx3,j) = 0;
% end
% 
% K(fx1,ma+1) = 1;
% K(fx2,ma+1) = 1;
% K(fx3,ma+1) = 1;
% R(fx1,1) = 0;
% R(fx2,1) = 0;
% R(fx3,1) = 0;

%-------------------------------------%
%   isostatic slab fixed at the top   %
%-------------------------------------%

% base = 2*(2*Tri(1).m+1);
% fx1 = 2*Tri(1).n*base + 1;
% fx2 = 2*Tri(1).n*base + 2;
% fx3 = (2*Tri(1).n + 1)*base;
% fxnodes(fx1,1) = 0;
% fxnodes(fx2,1) = 0;
% fxnodes(fx3,1) = 0;
% 
% for j=1:(2*ma+1)
%     K(fx1,j)=0;
%     K(fx2,j)=0;
%     K(fx3,j)=0;
% end
% 
% K(fx1,ma+1)=1;
% K(fx2,ma+1)=1;
% K(fx3,ma+1)=1;
% R(fx1,1)=0;
% R(fx2,1)=0;
% R(fx3,1)=0;

%-------------------------------------------------------%
%   isostatic beam fixed at the bottom and at the top   %
%-------------------------------------------------------%

% base = 2*(2*Tri(1).m+1);
% fx1 = 2*Tri(1).n*base + 1;
% fx2 = 2*Tri(1).n*base + 2;
% fx3 = 2;
% %fx3 = (2*Tri(1).n + 1)*base;
% fxnodes(fx1,1) = 0;
% fxnodes(fx2,1) = 0;
% fxnodes(fx3,1) = 0;
% 
% for j=1:(2*ma+1)
%     K(fx1,j)=0;
%     K(fx2,j)=0;
%     K(fx3,j)=0;
% end
% 
% K(fx1,ma+1)=1;
% K(fx2,ma+1)=1;
% K(fx3,ma+1)=1;
% R(fx1,1)=0;
% R(fx2,1)=0;
% R(fx3,1)=0;

%----------------------------------------%
%   isostatic slab fixed at the middle   %
%----------------------------------------%

base = 2*(2*Tri(1).m+1);
fx1 = Tri(1).n*base + 1;
fx2 = Tri(1).n*base + 2;
fx3 = (Tri(1).n + 1)*base;
fxnodes(fx1,1) = 0;
fxnodes(fx2,1) = 0;
fxnodes(fx3,1) = 0;

for j=1:(2*ma+1)
    K(fx1,j)=0;
    K(fx2,j)=0;
    K(fx3,j)=0;
end

K(fx1,ma+1)=1;
K(fx2,ma+1)=1;
K(fx3,ma+1)=1;
R(fx1,1)=0;
R(fx2,1)=0;
R(fx3,1)=0;

end