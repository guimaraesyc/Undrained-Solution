% CalcStiff calculates the updated stiffness matrix (K)

function K = CalcStiff_NL(delU1,Tri,dN,count)

K = zeros(12);                      % tangent stiffness matrix
K_NL = zeros(12);                   % linear and non linear part of stiffness matrix
K_sigma = zeros(12);                % geometric part of stiffness matrix
BL_real = zeros(3,12);              % linear part of kinematic matrix (real displacement)
BL_virtual = zeros(12,3);           % linear part of kinematic matrix (virtual displacement)
G = zeros(4,12);                    % derivative shape function matrix
Atheta = zeros(3,4);
H = zeros(3,4);
HAtheta = zeros(3,4);
delE = zeros(3,12);

S2 = zeros(4);
F = zeros(2);                % deformation gradient
Cgreen = zeros(2);           % right Cauchy-Green deformation tensor
Grad_NX = zeros(2,1);        % the derivate of shape function  with respect to reference coordinates
GradNx = zeros(2,1);         % the derivate of shape function  with respect current coordinates


for lnode=1:6
    gnode = Tri(1).ConM(count,lnode);
    u(lnode) = delU1(2*gnode-1);
    v(lnode) = delU1(2*gnode);
end


for ig=1:12 % loop for each Cowper´s integration point
    %for ig=1:7 % loop for each point of Gauss´s quadrature
    
    %---------------------------------------------------------------------
    %---------------------------------------------------------------------\
    G = zeros(4,12);                    % derivative shape function matrix
    Atheta = zeros(3,4);
    H = zeros(3,4);
    HAtheta = zeros(3,4);
    delE = zeros(3,12);
    K_NL = zeros(12);
    K_sigma = zeros(12);
    F = zeros(2);
    %---------------------------------------------------------------------\
    %---------------------------------------------------------------------\
    
    
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
    
    % Jacobian matrix (global to local coordinates)
    
    Jdxdcsi = dxdcsi(csi,eta,Tri,dN,count);
    Jdxdeta = dxdeta(csi,eta,Tri,dN,count);
    Jdydcsi = dydcsi(csi,eta,Tri,dN,count);
    Jdydeta = dydeta(csi,eta,Tri,dN,count);
    
    % Determinant of Jacobian matrix (|J|)
    
    WJdet = Jdxdcsi*Jdydeta - Jdxdeta*Jdydcsi;
    
    %--------------------------------------------------------------
    % Inverse Jacobian matrix (local to global coordinates)
    % J-¹=(1 / det J)*Jbar
    % Jbar: adjugate/adjunct matrix (transponse of cofactor matrix)
    %--------------------------------------------------------------
    
    Jdcsidx = Jdydeta/WJdet;
    Jdetadx = -Jdydcsi/WJdet;
    Jdcsidy = -Jdxdeta/WJdet;
    Jdetady = Jdxdcsi/WJdet;
    
    % weight of integration points (0.5 to fix the triangle area)
    %W_Jdet = 0.5*Tri(1).w_gp(ig)*W_Jdet;
    
    %{
    --------------------------------------------------------
    linear strain matrix (matrix BL)
    epsilon = B_L*u
    ------------------------------------------------------
    {epsilon x}    [ delta/delta x        0       ]   {ui}
    {epsilon y} =  [      0          delta/delta y] * {vi}
    {gama xy  }    [ delta/delta y   delta/delta x]
    ------------------------------------------------------
    %}
    
    WJdet = 0.5*Tri(1).w_cowper(ig)*WJdet;
    WJX_det = 0.5*Tri(1).w_cowper(ig)*WJX_det;
    %WJdet = 0.5*Tri(1).w_gp(ig)*WJdet;
    %WJX_det = 0.5*Tri(1).w_gp(ig)*WJX_det;
    
    %     for i=1:12
    %         BL_virtual(i,1) = mod(i,2)*(Jdcsidx*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetadx*dN(floor(i/2+0.5)).deta(csi,eta));
    %         BL_virtual(i,2) = mod((i+1),2)*(Jdcsidy*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetady*dN(floor(i/2+0.5)).deta(csi,eta));
    %         BL_virtual(i,3) = mod(i,2)*(Jdcsidy*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetady*dN(floor(i/2+0.5)).deta(csi,eta)) + mod((i+1),2)*(Jdcsidx*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetadx*dN(floor(i/2+0.5)).deta(csi,eta));
    %         for j=1:12
    %             BL_real(1,j) = mod(j,2)*(Jdcsidx*dN(floor(j/2+0.5)).dcsi(csi,eta) + Jdetadx*dN(floor(j/2+0.5)).deta(csi,eta));
    %             BL_real(2,j) = mod((j+1),2)*(Jdcsidy*dN(floor(j/2+0.5)).dcsi(csi,eta) + Jdetady*dN(floor(j/2+0.5)).deta(csi,eta));
    %             BL_real(3,j) = mod(j,2)*(Jdcsidy*dN(floor(j/2+0.5)).dcsi(csi,eta) + Jdetady*dN(floor(j/2+0.5)).deta(csi,eta)) + mod((j+1),2)*(Jdcsidx*dN(floor(j/2+0.5)).dcsi(csi,eta) + Jdetadx*dN(floor(j/2+0.5)).deta(csi,eta));
    %             for k=1:3
    %                 for l=1:3
    %                     K_L(i,j) = K_L(i,j) + WJdet*(BL_virtual(i,k)*Tri(1).C(k,l)*BL_real(l,j));
    %                 end
    %             end
    %         end
    %     end
    
    
    %----------------%
    %   NONLINEAR    %
    %----------------%
    
    for i=1:12
        
        G(1,i) = mod(i,2)*(Jdcsi_dX*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdeta_dX*dN(floor(i/2+0.5)).deta(csi,eta));
        G(2,i) = mod(i,2)*(Jdcsi_dY*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdeta_dY*dN(floor(i/2+0.5)).deta(csi,eta));
        G(3,i) = mod((i+1),2)*(Jdcsi_dX*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdeta_dX*dN(floor(i/2+0.5)).deta(csi,eta));
        G(4,i) = mod((i+1),2)*(Jdcsi_dY*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdeta_dY*dN(floor(i/2+0.5)).deta(csi,eta));
        
        %         G(1,i) = mod(i,2)*(Jdcsidx*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetadx*dN(floor(i/2+0.5)).deta(csi,eta));
        %         G(2,i) = mod(i,2)*(Jdcsidy*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetady*dN(floor(i/2+0.5)).deta(csi,eta));
        %         G(3,i) = mod((i+1),2)*(Jdcsidx*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetadx*dN(floor(i/2+0.5)).deta(csi,eta));
        %         G(4,i) = mod((i+1),2)*(Jdcsidy*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetady*dN(floor(i/2+0.5)).deta(csi,eta));
    end
    
    H(1,1) = 1;
    H(2,4) = 1;
    H(3,2) = 1;
    H(3,3) = 1;
    
    Atheta(1,1) = Jdcsi_dX*dudcsi(u,csi,eta,Tri,dN,count) + Jdeta_dX*dudeta(u,csi,eta,Tri,dN,count);
    Atheta(3,2) = Atheta(1,1);
    Atheta(2,2) = Jdcsi_dY*dudcsi(u,csi,eta,Tri,dN,count) + Jdeta_dY*dudeta(u,csi,eta,Tri,dN,count);
    Atheta(3,1) = Atheta(2,2);
    Atheta(1,3) = Jdcsi_dX*dvdcsi(v,csi,eta,Tri,dN,count) + Jdeta_dX*dvdeta(v,csi,eta,Tri,dN,count);
    Atheta(3,4) = Atheta(1,3);
    Atheta(2,4) = Jdcsi_dY*dvdcsi(v,csi,eta,Tri,dN,count) + Jdeta_dY*dvdeta(v,csi,eta,Tri,dN,count);
    Atheta(3,3) = Atheta(2,4);
    
    %     for i=1:3
    %         for j=1:4
    %             HAtheta(i,j) = H(i,j) + Atheta(i,j);
    %         end
    %     end
    %
    %     for i=1:3
    %         for j=1:12
    %             for k=1:4
    %                 delE(i,j) = delE(i,j) + HAtheta(i,k)*G(k,j);
    %             end
    %         end
    %     end
    
    HAtheta = H + Atheta;
    delE = HAtheta*G;
    K_NL = (delE'*Tri(1).C*delE);
    
    %     for i=1:12
    %         for j=1:12
    %             for k=1:3
    %                 for l=1:3
    %                     K_NL(i,j) = K_NL(i,j) + (delE(k,i)*Tri(1).C(k,l)*delE(l,j));
    %                 end
    %             end
    %         end
    %     end
    
    
    %------------------------------%
    %   DEFORMATION GRADIENT (F)   %
    %------------------------------%
    
    for lnode=1:6
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
    Cgreen = F'*F;
    %     for i=1:2
    %         for j=1:2
    %             for k=1:2
    %                 Cgreen(i,j) = Cgreen(i,j) + F(k,i)*F(k,j);
    %             end
    %         end
    %     end
    Egreen = 0.5*(Cgreen-eye(2));                                   % Green-Lagrange strain tensor
    Tstress = Tri(1).C*[Egreen(1,1);Egreen(2,2);2*Egreen(1,2)];     % 2nd PK in Voigt notation
    Spk = [Tstress(1,1),Tstress(3,1);Tstress(3,1),Tstress(2,1)];
    
    S2 = zeros(4);
    S2(1,1) = Spk(1,1);
    S2(1,2) = Spk(1,2);
    S2(2,1) = Spk(2,1);
    S2(2,2) = Spk(2,2);
    S2(3,3) = Spk(1,1);
    S2(3,4) = Spk(1,2);
    S2(4,3) = Spk(2,1);
    S2(4,4) = Spk(2,2);
    
    
    K_sigma = K_sigma + G'*S2*G;
    K = K + WJX_det*(K_NL + K_sigma);
    
    %     for i=1:12
    %         for j=1:12
    %             for k=1:4
    %                 for l=1:4
    %                     K_sigma(i,j) = K_sigma(i,j) + G(k,i)*S2(k,l)*G(l,j);
    %                 end
    %             end
    %         end
    %     end
    %
    %     for i=1:12
    %         for j=1:12
    %             K(i,j) = K(i,j) + WJX_det*(K_NL(i,j) + K_sigma(i,j));
    %         end
    %     end
end
end