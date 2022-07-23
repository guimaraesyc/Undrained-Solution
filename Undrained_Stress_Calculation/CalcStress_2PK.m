% this function calculates and plots the equivalent nodal stresses

function S2PKG = CalcStress_2PK(Tri,dN)
format long

S2PK = zeros(Tri(1).nnodes(1),1);
S2PK_count = zeros(Tri(1).nnodes(1),1);
S2PKG = zeros(Tri(1).nel,7);

I = eye(2);

for count=1:Tri(1).nel
    
    poly_count = 1;
    A = zeros(6);
    b = zeros(6,1);
    C = zeros(6,1);
    
    for ig=1:7  % 7 Gauss Legendre quadrature points
        csi = Tri(1).csi_gp(ig);
        eta = Tri(1).eta_gp(ig);
        
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
        Grad_NX = zeros(2,1);        % the derivate of shape function  with respect to reference coordinates
        
        for lnode=1:6
            Grad_NX(1,1) = (Jdcsi_dX*dN(lnode).dcsi(csi,eta) + Jdeta_dX*dN(lnode).deta(csi,eta));
            Grad_NX(2,1) = (Jdcsi_dY*dN(lnode).dcsi(csi,eta) + Jdeta_dY*dN(lnode).deta(csi,eta));
            
            F(1,1) = F(1,1) + Tri(count).x(lnode)*Grad_NX(1,1);
            F(1,2) = F(1,2) + Tri(count).x(lnode)*Grad_NX(2,1);
            F(2,1) = F(2,1) + Tri(count).y(lnode)*Grad_NX(1,1);
            F(2,2) = F(2,2) + Tri(count).y(lnode)*Grad_NX(2,1);
        end
        
        %------------------------------------------------------%
        %   DEFORMATION GRADIENT DUE TO IMPOSED STRAIN(Fres)   %
        %------------------------------------------------------%
        
        Fres = zeros(2);
        temp = 0;
        for lnode=1:6
            temp = temp + Tri(count).X(lnode)*dN(lnode).N(csi,eta);
        end
        
        T = Fourier(temp);
        Fres = Tri(count).alfa*(T-20) + I;
        
        
        %---------------------------------------%
        %   2ND PIOLA KIRCHHOFF STRESS TENSOR   %
        %---------------------------------------%
        
        Egreen = 0.5*(F'*F-I);                                                    % Green-Lagrange strain tensor
        %Eres = Fres;
        Eres = 0.5*(Fres'*Fres-I);                                                % residual strain tensor
        %Emec = Egreen;
        Emec = Egreen - Eres;                                                     % mechanical strain tensor
        Tstress = Tri(1).C*[Emec(1,1);Emec(2,2);2*Emec(1,2)];                     % 2nd PK in Voigt notation
        S = [Tstress(1,1),Tstress(3,1);Tstress(3,1),Tstress(2,1)];                % 2nd PK in matricial notation
        detF = F(1,1)*F(2,2) - F(1,2)*F(2,1);
        sigma = (1/detF)*F*S*F';                                                  % Cauchy stress
        S2PKG(count,ig) = sigma(2,2);
%         F_3D = [F(1,1),F(1,2),0;F(2,1),F(2,2),0;0,0,1];
%         temp = Tri(1).C(1,2)*(Emec(1,1)+Emec(2,2));
%         S_3D = [S(1,1),S(1,2),0;S(2,1),S(2,2),0;0,0,temp];
%         sigma_3D = (1/detF)*F_3D*S_3D*F_3D';                                                  % Cauchy stress
        
        % Ax=B, a(csi)2 + b(csi) + c(csi*eta) + d(eta) + e(eta)2 + f = S
        polynomial = zeros(1,6);
        if ig ~= 1
            polynomial = [(csi*csi),csi,(csi*eta),eta,(eta*eta),1];
            A(poly_count,:) = polynomial;
            %C(poly_count,1) = S(2,1);                                   % choose the component of stress tensor
            C(poly_count,1) = sigma(2,2);
            %C(poly_count,1) = sigma_3D(2,1);
            poly_count = poly_count + 1;
        end
    end
    
    % solve the linear system Ab = C
    b = gaussel(A,C);
    
    % calculate the stress contribution at each node in global coordinates
    for lnode=1:6
        csi = Tri(1).csi_nod(1,lnode);
        eta = Tri(1).eta_nod(1,lnode);
        gnode = Tri(1).ConM(count,lnode);
        S2PK(gnode,1) = S2PK(gnode,1) + b(1,1)*csi*csi + b(2,1)*csi + b(3,1)*csi*eta + b(4,1)*eta + b(5,1)*eta*eta + b(6,1);
        S2PK_count(gnode,1) = S2PK_count(gnode,1) + 1;
    end
end
% mean nodal stresses
for i=1:Tri(1).nnodes(1)
    S2PK(i,1) = S2PK(i,1) / S2PK_count(i,1);
end

% plot the results
gnode=0;
for i=1:(2*Tri(1).n+1)
    for j=1:(2*Tri(1).m+1)
        gnode = j+(i-1)*(2*Tri(1).m+1);
        %xx(i,j) = Tri(1).gnodes_0(2*gnode-1);
        %yy(i,j) = Tri(1).gnodes_0(2*gnode);
        xx(i,j) = Tri(1).gnodes(2*gnode-1);
        yy(i,j) = Tri(1).gnodes(2*gnode);
        zz(i,j) = S2PK(gnode,1);
    end
end
%contour(xx,yy,zz);
%h=contour(xx,yy,zz);
h=pcolor(xx,yy,zz);
set(h,'EdgeColor','none');
shading interp    % interpolate colors across lines and faces 
legend;
colorbar;
%lim = caxis;
%caxis([-5*10^6 5*10^6])
%axis([0, 0.4, 0, 0.4]);
%caxis([-10^5 10^5])
axis([-0.03, 0.06, 0, 1.03]);
end