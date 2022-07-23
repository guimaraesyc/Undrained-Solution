function K = Calcstiff_WI(Tri,dN,count)

K = zeros(12);

for gauss_point=1:7
    
csi = Tri(1).csi_gp(gauss_point);
eta = Tri(1).eta_gp(gauss_point);

JdX_dcsi = 0;
JdY_dcsi = 0;
JdX_deta = 0;
JdY_deta = 0;
for i=1:6
JdX_dcsi = JdX_dcsi + Tri(count).X(i)*dN(i).dcsi(csi,eta);
JdY_dcsi = JdY_dcsi + Tri(count).Y(i)*dN(i).dcsi(csi,eta);
JdX_deta = JdX_deta + Tri(count).X(i)*dN(i).deta(csi,eta);
JdY_deta = JdY_deta + Tri(count).Y(i)*dN(i).deta(csi,eta);
end
det_J = JdX_dcsi*JdY_deta - JdY_dcsi*JdX_deta;
Jdcsi_dX = JdY_deta/det_J;
Jdeta_dX = -JdY_dcsi/det_J;
Jdcsi_dY = -JdX_deta/det_J;
Jdeta_dY = JdX_dcsi/det_J;

% real part
B_real = zeros(3,12);
for i=1:12
temp = ceil(i/2);
B_real(1,i) = mod(i,2)*(Jdcsi_dX*dN(temp).dcsi(csi,eta) + Jdeta_dX*dN(temp).deta(csi,eta));
B_real(2,i) = mod(i+1,2)*(Jdcsi_dY*dN(temp).dcsi(csi,eta) + Jdeta_dY*dN(temp).deta(csi,eta));
B_real(3,i) = mod(i,2)*(Jdcsi_dY*dN(temp).dcsi(csi,eta) + Jdeta_dY*dN(temp).deta(csi,eta)) + mod(i+1,2)*(Jdcsi_dX*dN(temp).dcsi(csi,eta) + Jdeta_dX*dN(temp).deta(csi,eta));
end

% real virtual
B_virtual = zeros(12,3);
for i=1:12
for j=1:3
B_virtual(i,j) = B_real(j,i);
end
end

% Constitutive matrix
A = Tri(1).E/((1-2*Tri(1).ni)*(1+Tri(1).ni));
Tri(1).C = zeros(3);
Tri(1).C(1,1) = A*(1-Tri(1).ni);
Tri(1).C(1,2) = A*Tri(1).ni;
Tri(1).C(1,3) = 0;
Tri(1).C(2,1) = A*Tri(1).ni;
Tri(1).C(2,2) = A*(1-Tri(1).ni);
Tri(1).C(2,3) = 0;
Tri(1).C(3,1) = 0;
Tri(1).C(3,2) = 0;
Tri(1).C(3,3) = A*(0.5-Tri(1).ni);

WJ_det = 0.5*Tri(1).w_gp(gauss_point)*det_J;
for i=1:12
    for j=1:12
        for k=1:3
            for l=1:3
            K(i,j) = K(i,j) + WJ_det*B_virtual(i,k)*Tri(1).C(k,l)*B_real(l,j);
            end
        end
    end
end
end
end