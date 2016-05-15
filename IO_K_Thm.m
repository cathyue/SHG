% This is to calculate the nonlinear coupled mode equations involving Kerr
% Thermal nonlinearities
para;

% tunable:
dw = 0;
s0 = sqrt(1e-3);    %input 1mW
ke = [3e6; 6e6];

% material para:
ko = [3e6; 6e6];%Hz, corresponding to Q~2e8
Q = 2*pi*c0./(lam0.*(ko+ke));

func_final = 