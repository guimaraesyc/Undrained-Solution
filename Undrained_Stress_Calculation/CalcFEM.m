
% Displacement based Finite Element Method

function [delU,Tri] = CalcFEM (Tri,dN,Fe)
format long


%-----------------%
%   DECLARATION   %
%-----------------%

Fi = zeros(length(Fe),1);        % internal nodal forces (global nodes)
Residual = zeros(length(Fe),1);  % out of balance vector in Newton-Raphson iteraction process
delU = zeros(length(Fe),1);      % displacement due to applied load (global nodes)
delU1 = zeros(length(Fe),1);

Residual = Fe - Fi;
t=0;
error = 1;


%------------------------------%
%   NEWTON-RAPHSON PROCEDURE   %
%------------------------------%

for erro=1:1
    %while (error > 10^-15)
    t=t+1
    Fi = zeros(length(Fe),1);
    
    Kgel = Assmbl_BM(delU1,Tri,dN);
    [Kgel,Residual,fxnodes] = CalcBC_BM(Tri,Kgel,Residual);
    
    delU = Gaussel_BM (Kgel,Residual);
    delU1 = delU1 + delU;
    
    sqrt(delU'*delU)
    % update the currrent coordinates:  x,y (t+1) = x,y (t) + delta U,V
    factor = 1;
    Tri(1).gnodes = Tri(1).gnodes + factor*delU;
    for count=1:Tri(1).nel
        for lnode=1:6
            gnode = Tri(1).ConM(count,lnode);
            Tri(count).x(lnode) = Tri(1).gnodes(2*gnode-1,1);
            Tri(count).y(lnode) = Tri(1).gnodes(2*gnode,1);
        end
    end
    
    Fi = CalcFint(Tri,dN);
    Residual = Fe - Fi;
    t;
    
    % This part evaluate the error at each increment in Newton-Raphson
    % iteration
    error = 0;
    tol1=0;
    tol2=0;
    for count=1:length(Fe)
        tol1 = tol1 + fxnodes(count,1)*sqrt(Residual(count,1)*Residual(count,1));
        tol2 = tol2 + fxnodes(count,1)*sqrt(Fi(count,1)*Fi(count,1));
    end
    error = tol1/tol2
end

% total displacement vector (u = x - X)
for count=1:Tri(1).nel(1)
    for lnode=1:6
        gnode = Tri(1).ConM(count,lnode);
        delU(2*gnode-1) = Tri(count).x(lnode) - Tri(count).X(lnode);
        delU(2*gnode) = Tri(count).y(lnode) - Tri(count).Y(lnode);
    end
end

%--------------------------%
%   PLOT TRIANGULAR MESH   %
%--------------------------%

% for count=1:Tri(1).nel(1)
%     for j=1:3
%         temp = Tri(1).ConM(count,j);
%         T1(count,j) = temp;
%         xx(temp) = Tri(1).gnodes(2*temp-1);
%         yy(temp) = Tri(1).gnodes(2*temp);
%     end
% end
% triplot(T1,xx,yy)
% axis([-0.35, 0.35, -0.35, 0.35])
%axis([-0.05, 0.35, -0.05, 0.35])

% plot the results in a color graphic
%[xx,yy,zz] = Pprod(Tri)
end