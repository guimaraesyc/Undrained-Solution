
function temp = dvdeta (v,csi,eta,Tri,dN,count)
format long
	
temp = 0;
for i=1:6
temp = temp + dN(i).deta(csi,eta)*v(i);
end

end
