% This function creates the database used in FEM analysis

function [dN,Tri] = CalcDB (Tri)
format long

% dN.N: shape functions
% dN.dcsi, dN.deta: derivative of shape functions with respect to local
% coordinates
% Tri.eta_gp,csi_gp,w_gp: Gauss integration points and weights
% Tri.alfa: calcite thermal expansion coef
% Tri.C: elastic matrix (E,ni)
% Tri.(X,Y): coordinates of initial configuration

%----------------------------------%
%   SHAPE FUNCTION & DERIVATIVES   %
%----------------------------------%

% Shape functions
dN(1).N = @(csi,eta) ((1.-(csi+eta))*(1.-2.*(csi+eta)));
dN(2).N = @(csi,eta) (csi*(2*csi-1.)+0*eta);
dN(3).N = @(csi,eta) (eta*(2*eta-1.)+0*csi);
dN(4).N = @(csi,eta) (4.*csi*(1.-(csi+eta)));
dN(5).N = @(csi,eta) (4.*eta*(1.-(csi+eta)));
dN(6).N = @(csi,eta) (4.*eta*csi);
% Derivative of shape function in relation to csi
dN(1).dcsi = @(csi,eta) (4*(csi+eta)-3);       % d_W_0_d_csi
dN(2).dcsi = @(csi,eta) (4*csi-1 + 0*eta);     % d_W_1_d_csi
dN(3).dcsi = @(csi,eta) (0*csi + 0*eta);       % d_W_2_d_csi
dN(4).dcsi = @(csi,eta) (4*(1-2*csi-eta));     % d_W_3_d_csi
dN(5).dcsi = @(csi,eta) (0*csi - 4*eta);       % d_W_4_d_csi
dN(6).dcsi = @(csi,eta) (0*csi + 4*eta);       % d_W_5_d_csi
% Derivative of shape function in relation to eta
dN(1).deta = @(csi,eta) (4*(csi+eta)-3);       % d_W_0_d_eta
dN(2).deta = @(csi,eta) (0*csi + 0*eta);       % d_W_1_d_eta
dN(3).deta = @(csi,eta) (0*csi + 4*eta-1);     % d_W_2_d_eta
dN(4).deta = @(csi,eta) (-4*csi + 0*eta);      % d_W_3_d_eta
dN(5).deta = @(csi,eta) (4*(1-(csi+2*eta)));   % d_W_4_d_eta
dN(6).deta = @(csi,eta) (4*csi + 0*eta);       % d_W_5_d_eta


%----------------------------------------%
%   GAUSS QUADRATURE - QUINTIC SPLINES   %
%----------------------------------------%

alpha_1 = 0.059715871789768;
alpha_2 = 0.797426985353087;
beta_1 = 0.470142064105115;
beta_2 = 0.101286507323456;
wg_1 = 0.225;
wg_2 = 0.132394152788506;
wg_3 = 0.125939180544827;
Tri(1).csi_nod = [0, 1, 0, 0.5, 0, 0.5];
Tri(1).eta_nod = [0, 0, 1, 0, 0.5, 0.5];
Tri(1).csi_gp = [1/3, alpha_1, beta_1, beta_1, alpha_2, beta_2, beta_2];
Tri(1).eta_gp = [1/3, beta_1, alpha_1, beta_1, beta_2, alpha_2, beta_2];
Tri(1).w_gp = [wg_1, wg_2, wg_2, wg_2, wg_3, wg_3, wg_3];


%--------------------------------%
%   COOPER FORMULA - 12 points   %
%--------------------------------%

alpha_1 = 0.873821971016996;
beta_1 = 0.063089014491502;
gama_1 = 0.063089014491502;
wg_1 = 0.050844906370207;

alpha_2 = 0.501426509658179;
beta_2 = 0.249286745170910;
gama_2 = 0.249286745170911;
wg_2 = 0.116786275726379;

alpha_3 = 0.636502499121399;
beta_3 = 0.310352451033785;
gama_3 = 0.053145049844816;
wg_3 = 0.082851075618374;

Tri(1).csi_cowper = [alpha_1, beta_1, beta_1, alpha_2, beta_2, beta_2, alpha_3, alpha_3, beta_3, beta_3, gama_3, gama_3];
Tri(1).eta_cowper = [beta_1, alpha_1, beta_1, beta_2, alpha_2, beta_2, beta_3, gama_3, alpha_3, gama_3, alpha_3, beta_3];
Tri(1).w_cowper = [wg_1, wg_1, wg_1, wg_2, wg_2, wg_2, wg_3, wg_3, wg_3, wg_3, wg_3, wg_3];


%-------------------------%
%   TRIANGULAR ELEMENTS   %
%-------------------------%


for count=1:Tri(1).nel
    
    % Elastic parameters
    Tri(count).E = 52.4*10^9;                   % Pa  (Ferrero et al, 2009)
    Tri(count).ni = 0.16;                       % (Ferrero et al, 2009)
    %Tri(count).E = 1;
    %Tri(count).ni = 0.25;
    
    % Thermal expansion coefficient
    Tri(count).alfa = zeros(2);
    Tri(count).alfa(1,1) = 5.9*10^-6;    % °C-1  (Ferrero et al, 2009)
    Tri(count).alfa(2,2) = 5.9*10^-6;
    
    % Triangular element
    Tri(count).X = zeros(6,1);
    Tri(count).Y = zeros(6,1);
    Tri(count).node = zeros(6,1);
    Tri(count).x = zeros(6,1);
    Tri(count).y = zeros(6,1);
    
    % initial displacement
    Tri(count).u = zeros(6,1);
    Tri(count).v = zeros(6,1);
    
    % Elastic matrix
    Tri(count).C=zeros(3);
end

for count=1:Tri(1).nel
    Tri(count).C = CalcC(Tri); % Elastic matrix
    for lnode=1:6
        temp = Tri(1).ConM(count,lnode);
        Tri(count).X(lnode) = Tri(1).gnodes_0((2*temp-1));
        Tri(count).Y(lnode) = Tri(1).gnodes_0((2*temp));
        Tri(count).x(lnode) = Tri(1).gnodes_0((2*temp-1));
        Tri(count).y(lnode) = Tri(1).gnodes_0((2*temp));
    end
end

end
