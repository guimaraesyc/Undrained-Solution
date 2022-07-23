function temp = dydcsi (csi,eta,Tri,dN,count)
format long
	
temp = 0;
for i=1:6
temp = temp + dN(i).dcsi(csi,eta)*Tri(count).y(i);
end
end