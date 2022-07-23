function temp = dY_dcsi(csi,eta,Tri,dN,count)
format long

temp = 0;
for i=1:6
temp = temp + dN(i).dcsi(csi,eta)*Tri(count).Y(i);
end
end