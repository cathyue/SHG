load sample.mat;
% tunable:
dw = 0;
s0 = sqrt(1e-3);    %input 1mW
% lam1 = 1.55e-6;     %m
% lam2 = lam1/2;  %suppose on resonance
lam0 = [lam10; lam20];
n0 = [n10; n20];
w0 = [w10; w20];
l0 = [l0; 2*l0];

% input para:
delt = 7.0649e+09; %from detuning.m
kmm = -2.0953e-10 - 8.8476e+05i; %from kappa_cal.m
% n1 = n10;   %from n_lam.m
% n2 = n20;   

% constant:
c0 = 299792458; %m/s
Z0 = 376.73;    %free space resistance
epsi0 = 8.854188e-12;   %F/m

% material para:
% ko1 = 3e6; %Hz
% ke1 = 3e6; %Hz
% ko2 = 6e6; %Hz
% ke2 = 6e6; 
ko = [3e6; 6e6];%Hz, corresponding to Q~2e8
ke = [3e6; 6e6];
Q = 2*pi*c0./(lam0.*(ko+ke));
% Q2 = 2*pi*c0/(lam20*(ko2+ke2));

n2 = [2.79e-20; 2.48e-20]; %m2/W, Review and assessment of measured values of the nonlinear refractive-index coefficient of fused silica David Milam
kai3 = n2.*4.*n0.^2.*epsi0*c0/3;   %m2/V2
% kai32 = n2_2*4*n20^2*epsi0*c0/3;    %m2/V2

dndT = 8.7e-6;  %1/K
rho = 2200; %kg/m3
C = 740;    %J/(kgK)
D = 9.5e-7; %m2/s

Qab = [7e8; 2e10];  %from Rokhsari_APL_2004, Fig.2

b = [1.7e-6;sqrt(1.08^2+1)*1e-6];    %m, estimated from phil(l,x)

dthet = 2*D./(b.^2);







