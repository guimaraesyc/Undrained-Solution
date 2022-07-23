% This function calculates the nodal forces due to an imposed displacements

function Fi = CalcFint(Tri,dN)
format long

Fi = zeros(2*Tri(1).nnodes,1);    % internal nodal forces (local nodes)

%-------------------------%
%   NODAL FORCES (Fint)   %
%-------------------------%

for count=1:Tri(1).nel
    for ig=1:12                       % 12 Cowper's integration points
    %for ig=1:7                       % 7 Gauss Legendre quadrature points
        
        csi = Tri(1).csi_cowper(ig);
        eta = Tri(1).eta_cowper(ig);
        %csi = Tri(1).csi_gp(ig);
        %eta = Tri(1).eta_gp(ig);
        
        % Jacobian matrix from local to reference coordinates
        JdX_dcsi = dX_dcsi(csi,eta,Tri,dN,count);
        JdY_dcsi = dY_dcsi(csi,eta,Tri,dN,count);
        JdX_deta = dX_deta(csi,eta,Tri,dN,count);
        JdY_deta = dY_deta(csi,eta,Tri,dN,count);
        
        WJX_det  = JdX_dcsi*JdY_deta - JdY_dcsi*JdX_deta;
        
        Jdcsi_dX =  JdY_deta/WJX_det;
        Jdcsi_dY = -JdX_deta/WJX_det;
        Jdeta_dX = -JdY_dcsi/WJX_det;
        Jdeta_dY =  JdX_dcsi/WJX_det;
        
        
        %------------------------------%
        %   DEFORMATION GRADIENT (F)   %
        %------------------------------%
        
        F = zeros(2);
        for lnode=1:6
            Grad_NX = zeros(2,1);        % the derivate of shape function  with respect to reference coordinates
            Grad_NX(1,1) = (Jdcsi_dX*dN(lnode).dcsi(csi,eta) + Jdeta_dX*dN(lnode).deta(csi,eta));
            Grad_NX(2,1) = (Jdcsi_dY*dN(lnode).dcsi(csi,eta) + Jdeta_dY*dN(lnode).deta(csi,eta));
            
            F(1,1) = F(1,1) + Tri(count).x(lnode)*Grad_NX(1,1);
            F(1,2) = F(1,2) + Tri(count).x(lnode)*Grad_NX(2,1);
            F(2,1) = F(2,1) + Tri(count).y(lnode)*Grad_NX(1,1);
            F(2,2) = F(2,2) + Tri(count).y(lnode)*Grad_NX(2,1);
        end
        
        
        %---------------------------------------%
        %   2ND PIOLA KIRCHHOFF STRESS TENSOR   %
        %---------------------------------------%
        
        Cgreen = F'*F;                                                  % right Cauchy-Green tensor
        Egreen = 0.5*(Cgreen-eye(2));                                   % Green-Lagrange strain tensor
        Tstress = Tri(1).C*[Egreen(1,1);Egreen(2,2);2*Egreen(1,2)];     % 2nd PK in Voigt notation
        Spk = [Tstress(1),Tstress(3);Tstress(3),Tstress(2)];            % 2nd PK in matrix notation
        Spk = F*Spk;
        
        
        %------------------%
        %   NODAL FORCES   %
        %------------------%
        
        W_Jdet  = 0.5*Tri(1).w_cowper(ig)*WJX_det;  % integrating over the initial volume (X,Y)
        %W_Jdet  = 0.5*Tri(1).w_gp(ig)*WJX_det;  % integrating over the initial volume (X,Y)
        
        for lnode=1:6
            Grad_NX = zeros(2,1);        % the derivate of shape function  with respect to reference coordinates
            Grad_NX(1,1) = (Jdcsi_dX*dN(lnode).dcsi(csi,eta) + Jdeta_dX*dN(lnode).deta(csi,eta));
            Grad_NX(2,1) = (Jdcsi_dY*dN(lnode).dcsi(csi,eta) + Jdeta_dY*dN(lnode).deta(csi,eta));
            
            gnode = Tri(1).ConM(count,lnode);
            Fi((2*gnode-1),1) = Fi((2*gnode-1),1) + W_Jdet*(Spk(1,1)*Grad_NX(1,1) + Spk(1,2)*Grad_NX(2,1));
            Fi(2*gnode,1) = Fi(2*gnode,1) + W_Jdet*(Spk(2,1)*Grad_NX(1,1)+ Spk(2,2)*Grad_NX(2,1));
        end
    end
end
end