% Calculate the elastic matrix (2D - Plane strain state)

function C = CalcC (Tri)
format long
	
	C=zeros(3);
	C(1,1) = Tri(1).E*(1-Tri(1).ni)/((1+Tri(1).ni)*(1-2*Tri(1).ni));
	C(1,2) = Tri(1).E*Tri(1).ni/((1+Tri(1).ni)*(1-2*Tri(1).ni));
	C(1,3) = 0;
	C(2,1) = Tri(1).E*Tri(1).ni/((1+Tri(1).ni)*(1-2*Tri(1).ni));
	C(2,2) = Tri(1).E*(1-Tri(1).ni)/((1+Tri(1).ni)*(1-2*Tri(1).ni));
	C(2,3) = 0;
	C(3,1) = 0;
	C(3,2) = 0;
	C(3,3) = Tri(1).E/(2*(1+Tri(1).ni));
	
end


%{
{sigma x }      [C11    C12    0 ]       {epsilon x}
{sigma y }   =  [C21    C22    0 ]   *   {epsilon y}
{ tal xy }      [  0     0    C33]       { gama xy }
%}
