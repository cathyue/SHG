% IO_K_Thm.m
% This is to calculate the nonlinear coupled mode equations involving Kerr
% Thermal nonlinearities
para;

% tunable:
dw = 0;
s0 = sqrt(3.46e-4);    %input sqrt(W)
ke = [3e6; 6e6].*(2*pi);    %rad/s

% material para:
ko = [3e6; 6e6].*(2*pi);%rad/s, corresponding to Q~2e8
Q = 2*pi*c0./(lam0.*(ko+ke)./(2*pi));

a0 = sqrt(ke).*s0./((ko+ke)./2);
a0 = a0.*[1; 0];
a0(2) = kmm*a0(1)^2/((ko(2)+ke(2))/2);
% kmm = 0;
sweep = 3.5e10:2e8:4e10;
% B = 0.*B;
% delt = 0;
theta0 = [0; 0];
a0 = a0.*[1*(cos(theta0(1))+1i*sin(theta0(1))); 0.5*(cos(theta0(2))+1i*sin(theta0(2)))];
x = [real(a0(1)), imag(a0(1)), real(a0(2)), imag(a0(2))];
P1 = zeros(length(sweep),1);
P2 = zeros(length(sweep),1);
flag = zeros(length(sweep), 1);
options = optimoptions('fsolve','MaxFunEval', 1600,'MaxIter', 1600);
for kw = 1:length(sweep)
    dw = sweep(kw);
    %     if kw<84
    %         x0 = 0.*x;
    %     else
    %         x0 = x;
    %     end
    x0 = x;
    %     exitflag = -1;
    [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
    flag(kw) = exitflag;
    if exitflag <0
        x0 = 0.*x;
        [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
    end
    a = [func_final(1)+1i*func_final(2); func_final(3)+1i*func_final(4)];
    Pout = abs(a.*sqrt(ke)).^2;
    P1(kw) = Pout(1);
    P2(kw) = Pout(2);
end
% figure; plot(sweep, P1);