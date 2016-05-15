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

a0 = sqrt(ke).*s0./((ko+ke)./2);
a0 = a0.*[1;0];
kmm = 0;
% B = 0.*B;
theta0 = [0; 0];
a0 = a0.*[1*(cos(theta0(1))+1i*sin(theta0(1))); 0.5*(cos(theta0(2))+1i*sin(theta0(2)))];
x0 = [real(a0(1)), imag(a0(1)), real(a0(2)), imag(a0(2))];
func_final = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0);
a = [func_final(1)+1i*func_final(2); func_final(3)+1i*func_final(4)];
Pout = abs(a.*sqrt(ke)).^2;