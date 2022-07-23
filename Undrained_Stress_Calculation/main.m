clear;
clc;
format long


%--------------%
%   INC DATA   %
%--------------%
m = 20;                                           % row subdivision
n = 5;                                           % column subdivision

[Tri] = Rdata(m,n);

[dN,Tri] = CalcDB (Tri);
Tri(1).m = m;                                    % row subdivision
Tri(1).n = n;                                    % column subdivision

%------------------%
%   THERMAL LOAD   %
%------------------%

Fthermal = CalcFthermal(Tri,dN);

%------------%
%   SOLVER   %
%------------%

[delU,Tri] = CalcFEM(Tri,dN,Fthermal);

%-------------%
%   RESULTS   %
%-------------%

S2PKG = CalcStress_2PK(Tri,dN);