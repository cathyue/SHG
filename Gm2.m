function [ Gmm ] = Gm2( l, x, n )
%UNTITLED3 此处显示有关此函数的摘要
%   my notebook

Gmm = -n.^2.*jl(l-2, n.*x)...
    +jl(l-1, n.*x).*(n.*(2.*l+1)./x-(l+1).*n.^3/(x)) ...
    -jl(l, n.*x).*(l+1)./x.^2.*(1-n.^2).*(l+2) ...
    +n.^2/(hl(l, x)).*(hl(l-2, x).*jl(l, n.*x)+n.*hl(l-1, x).*jl(l-1, n.*x)...
                        - hl(l-1, x).*jl(l, n.*x).*l./x) ...
    -n.^2.*hl(l-1, x).^2.*jl(l, n.*x)./(hl(l, x).^2);


end

