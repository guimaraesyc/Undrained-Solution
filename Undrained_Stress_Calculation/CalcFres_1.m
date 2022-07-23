% This function calculates the nodal forces due to an imposed displacements

function Fres1 = CalcFres_1(Tri,dN)
format long

% imposed strain due to bowing
Eres1 = 0.02;       % imposed strain (internal face)
Eres2 = 0.03;       % imposed strain (external face)
I = eye(2);
E = zeros(2);
Fres1 = zeros(2*Tri(1).nnodes(1),1);

F1 = eye(2)*(1+Eres1);
F2 = eye(2)*(1+Eres2);

% E = [H + A_theta]*theta
% linear part of Green strain in matricial notation
H = zeros(3,4);
H(1,1) = 1;
H(2,4) = 1;
H(3,2) = 1;
H(3,3) = 1;

for count=1:Tri(1).nel
    
    for ig=1:7       % 7 Gauss Legendre quadrature points
        
        csi = Tri(1).csi_gp(ig);
        eta = Tri(1).eta_gp(ig);
        
        %------------------------------%
        %   DEFORMATION GRADIENT (F)   %
        %------------------------------%
        
        F = zeros(2);
        temp = 0;
        for i=1:6
            temp = temp + Tri(count).X(i)*dN(i).N(csi,eta);
        end
        % linear interpolation y = a + bx
        F(1,1) = F1(1,1) + (temp/Tri(1).b(1))*(F2(1,1)-F1(1,1));
        F(2,2) = F(1,1);
        
        % displacement gradient (D = F - I)
        D = zeros(2);
        D = F - eye(2);
        % displacement gradient in vectorized form (theta)
        theta = [D(1,1);D(1,2);D(2,1);D(2,2)];
        % nonlinear part
        A_theta = zeros(3,4);
        temp = [2,1,4,3];
        for j=1:4
            A_theta(1,j) = mod(j,2)*theta(j,1);
            A_theta(2,j) = mod((j+1),2)*theta(j,1);
            A_theta(3,j) = theta(temp(1,j),1);
        end
        
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
        
        G = zeros(4,12);
        for j=1:12
            temp = floor ((j+1)/2);
            G(1,j) = mod(j,2)*(Jdcsi_dX*dN(temp).dcsi(csi,eta) + Jdeta_dX*dN(temp).deta(csi,eta));
            G(2,j) = mod(j,2)*(Jdcsi_dY*dN(temp).dcsi(csi,eta) + Jdeta_dY*dN(temp).deta(csi,eta));
            G(3,j) = mod((j+1),2)*(Jdcsi_dX*dN(temp).dcsi(csi,eta) + Jdeta_dX*dN(temp).deta(csi,eta));
            G(4,j) = mod((j+1),2)*(Jdcsi_dY*dN(temp).dcsi(csi,eta) + Jdeta_dY*dN(temp).deta(csi,eta));
        end
        
        
        %---------------------------------------%
        %   2ND PIOLA KIRCHHOFF STRESS TENSOR   %
        %---------------------------------------%
        
        Cgreen = F'*F;                                                  % right Cauchy-Green tensor
        Egreen = 0.5*(Cgreen-eye(2));                                   % Green-Lagrange strain tensor
        S2pk = Tri(1).C*[Egreen(1,1);Egreen(2,2);2*Egreen(1,2)];        % 2nd PK in Voigt notation
        
        
        %---------------------------------------%
        %   2ND PIOLA KIRCHHOFF STRESS TENSOR   %
        %---------------------------------------%
        
        W_Jdet  = 0.5*Tri(1).w_gp(ig)*WJX_det;  % integrating over the initial volume (X,Y)
        del_E = (H + A_theta)*G*W_Jdet;                                        % increment of Green strain
        
        
        %---------------------------%
        %   NODAL FORCES (Fres,i)   %
        %---------------------------%
        temp = zeros(12,1);
        temp = del_E'*S2pk;
        
        for lnode=1:6
            % equivalent nodal forces in global coordinates
            gnode = Tri(1).ConM(count,lnode);
            Fres1((2*gnode-1),1) = Fres1((2*gnode-1),1) + temp(2*lnode-1,1);
            Fres1(2*gnode,1) = Fres1(2*gnode,1) + temp(2*lnode,1);
        end
    end
end
end