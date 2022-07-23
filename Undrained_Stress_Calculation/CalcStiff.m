% CalcStiff calculates the updated stiffness matrix (K)

function K = CalcStiff(Tri,dN,count)

K = zeros(12);                      % stiffness matrix
BL_real = zeros(3,12);              % linear part of kinematic matrix (real displacement)
BL_virtual = zeros(12,3);           % linear part of kinematic matrix (virtual displacement)

for ig=1:7 % loop for each point of Gauss´s quadrature
    
    csi = Tri(1).csi_gp(ig);
    eta = Tri(1).eta_gp(ig);
    
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
    
    WJdet = 0.5*Tri(1).w_gp(ig)*WJdet;
        for i=1:12
            BL_virtual(i,1) = mod(i,2)*(Jdcsidx*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetadx*dN(floor(i/2+0.5)).deta(csi,eta));
            BL_virtual(i,2) = mod((i+1),2)*(Jdcsidy*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetady*dN(floor(i/2+0.5)).deta(csi,eta));
            BL_virtual(i,3) = mod(i,2)*(Jdcsidy*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetady*dN(floor(i/2+0.5)).deta(csi,eta)) + mod((i+1),2)*(Jdcsidx*dN(floor(i/2+0.5)).dcsi(csi,eta) + Jdetadx*dN(floor(i/2+0.5)).deta(csi,eta));
            for j=1:12
                BL_real(1,j) = mod(j,2)*(Jdcsidx*dN(floor(j/2+0.5)).dcsi(csi,eta) + Jdetadx*dN(floor(j/2+0.5)).deta(csi,eta));
                BL_real(2,j) = mod((j+1),2)*(Jdcsidy*dN(floor(j/2+0.5)).dcsi(csi,eta) + Jdetady*dN(floor(j/2+0.5)).deta(csi,eta));
                BL_real(3,j) = mod(j,2)*(Jdcsidy*dN(floor(j/2+0.5)).dcsi(csi,eta) + Jdetady*dN(floor(j/2+0.5)).deta(csi,eta)) + mod((j+1),2)*(Jdcsidx*dN(floor(j/2+0.5)).dcsi(csi,eta) + Jdetadx*dN(floor(j/2+0.5)).deta(csi,eta));
                for k=1:3
                    for l=1:3
                        K(i,j) = K(i,j) + WJdet*(BL_virtual(i,k)*Tri(1).C(k,l)*BL_real(l,j));
                    end
                end
            end
        end
    
    %-----PI-----PI-----PI-----PI-----PI-----PI-----
    %-----PI-----PI-----PI-----PI-----PI-----PI-----
%     for i=1:12
%         Strain_i(1) = (i%2?0.:1.)*d_W_d_csi[int(i/2)](csi,eta)*JdcsidX + (i%2?0.:1.)*d_W_d_eta[int(i/2)](csi,eta)*JdetadX; //strain due to displ. i
%         Strain_i(2) = (i%2?1.:0.)*d_W_d_csi[int(i/2)](csi,eta)*JdcsidY + (i%2?1.:0.)*d_W_d_eta[int(i/2)](csi,eta)*JdetadY;
%         Strain_i(3) = (i%2?0.:1.)*d_W_d_csi[int(i/2)](csi,eta)*JdcsidY + (i%2?0.:1.)*d_W_d_eta[int(i/2)](csi,eta)*JdetadY
%                     + (i%2?1.:0.)*d_W_d_csi[int(i/2)](csi,eta)*JdcsidX + (i%2?1.:0.)*d_W_d_eta[int(i/2)](csi,eta)*JdetadX;
%         for j=1:12
%             Strain_j(1) = (j%2?0.:1.)*d_W_d_csi[int(j/2)](csi,eta)*JdcsidX + (j%2?0.:1.)*d_W_d_eta[int(j/2)](csi,eta)*JdetadX; //strain due to virtual displ.j
%             Strain_j(2) = (j%2?1.:0.)*d_W_d_csi[int(j/2)](csi,eta)*JdcsidY + (j%2?1.:0.)*d_W_d_eta[int(j/2)](csi,eta)*JdetadY;
%             Strain_j(3) = (j%2?0.:1.)*d_W_d_csi[int(j/2)](csi,eta)*JdcsidY + (j%2?0.:1.)*d_W_d_eta[int(j/2)](csi,eta)*JdetadY
%             +(j%2?1.:0.)*d_W_d_csi[int(j/2)](csi,eta)*JdcsidX + (j%2?1.:0.)*d_W_d_eta[int(j/2)](csi,eta)*JdetadX;
%             for k=1:12
%                 for l=1:12
%                     K(i,j) = K(i,j) + W_Jdet*Strain_i(k)*Strain_j(l)*Tri(1).C(k,l);
%                 end
%             end
%         end
%     end
    %-----PI-----PI-----PI-----PI-----PI-----PI-----
    %-----PI-----PI-----PI-----PI-----PI-----PI-----
end
end