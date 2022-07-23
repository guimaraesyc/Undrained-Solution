% This function calculates the nodal forces due to an imposed displacements

function Force_res = CalcFthermal(Tri,dN)
format long

I = eye(2);
Force_res = zeros(2*Tri(1).nnodes,1);

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
        temp = 0;
        for i=1:6
            temp = temp + Tri(count).X(i)*dN(i).N(csi,eta);
        end
        
        T = Fourier(temp);
        F = Tri(count).alfa*(T-20) + eye(2);
        
        
        %---------------------------------------%
        %   2ND PIOLA KIRCHHOFF STRESS TENSOR   %
        %---------------------------------------%
        
        Cgreen = F'*F;                                                           % right Cauchy-Green tensor
        Egreen = 0.5*(Cgreen-eye(2));                                            % Green-Lagrange strain tensor
        Tstress = Tri(1).C*[Egreen(1,1);Egreen(2,2);2*Egreen(1,2)];              % 2nd PK in Voigt notation
        Spk1 = [Tstress(1,1),Tstress(3,1);Tstress(3,1),Tstress(2,1)];             % 2nd PK in matrix notation
        Spk = F*Spk1;
        
%         Cgreen = zeros(2);
%         for i=1:2
%             for j=1:2
%                 for k=1:2
%                     Cgreen(i,j) = Cgreen(i,j) + F(k,i)*F(k,j);
%                 end
%             end
%         end
%         
%         Egreen = zeros(2);
%         Egreen(1,1) = 0.5*(Cgreen(1,1) - 1);
%         Egreen(1,2) = 0.5*(Cgreen(1,2) - 0);
%         Egreen(2,1) = 0.5*(Cgreen(2,1) - 0);
%         Egreen(2,2) = 0.5*(Cgreen(2,2) - 1);
%         E_voigt = zeros(3,1);
%         E_voigt(1,1) = Egreen(1,1);
%         E_voigt(2,1) = Egreen(2,2);
%         E_voigt(3,1) = Egreen(1,2) + Egreen(2,1);
%         
%         Tstress = zeros(3,1);
%         for i=1:3
%             for j=1:1
%                 for k=1:3
%                     Tstress(i,j) = Tstress(i,j) + Tri(1).C(i,k)*E_voigt(k,j);
%                 end
%             end
%         end
%         
%         Spk1 = zeros(2);
%         Spk1(1,1) = Tstress(1,1);
%         Spk1(1,2) = Tstress(3,1);
%         Spk1(2,1) = Tstress(3,1);
%         Spk1(2,2) = Tstress(2,1);
%         
%         Spk = zeros(2,2);
%         for i=1:2
%             for j=1:2
%                 for k=1:2
%                     Spk(i,j) = Spk(i,j) + F(i,k)*Spk1(k,j);
%                 end
%             end
%         end
        
        
        %---------------------------%
        %   NODAL FORCES (Fres,i)   %
        %---------------------------%
        
        for lnode=1:6
            Grad_NX = zeros(2,1);        % the derivate of shape function  with respect to reference coordinates (dNi / dX)
            Grad_NX(1,1) = (Jdcsi_dX*dN(lnode).dcsi(csi,eta) + Jdeta_dX*dN(lnode).deta(csi,eta));
            Grad_NX(2,1) = (Jdcsi_dY*dN(lnode).dcsi(csi,eta) + Jdeta_dY*dN(lnode).deta(csi,eta));
            
            W_Jdet  = 0.5*Tri(1).w_cowper(ig)*WJX_det;  % integrating over the initial volume (X,Y)
            %W_Jdet  = 0.5*Tri(1).w_gp(ig)*WJX_det;  % integrating over the initial volume (X,Y)
            
            % equivalent nodal forces in global coordinates
            gnode = Tri(1).ConM(count,lnode);
            Force_res((2*gnode-1),1) = Force_res((2*gnode-1),1) + W_Jdet*(Spk(1,1)*Grad_NX(1,1) + Spk(1,2)*Grad_NX(2,1));
            Force_res(2*gnode,1) = Force_res(2*gnode,1) + W_Jdet*(Spk(2,1)*Grad_NX(1,1)+ Spk(2,2)*Grad_NX(2,1));
        end
    end
end
end