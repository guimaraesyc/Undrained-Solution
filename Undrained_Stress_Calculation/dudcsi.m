
function temp = dudcsi (u,csi,eta,Tri,dN,count)
format long
	
temp = 0;
for i=1:6
temp = temp + dN(i).dcsi(csi,eta)*u(i);
end

end