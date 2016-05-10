function [ G ] = Gm( l,x,n )
%UNTITLED4 Summary of this function goes here
%   Kozyreff PRA 2008, eq. (35)

% G = hl(l-2, x)-l/x*hl(l-1,x)-(hl(l-1,x)-(l+1)/x*hl(l,x))*jl(l-1,x)/(n*(jl(l,x)))...
%     -hl(l,x)/n*((jl(l-2,x)-l/x*jl(l-1,x))*jl(l,x)-(jl(l-1,x)-(l+1)/x*jl(l,x))*jl(l-1,x))/(jl(l-1,x)^2)...
%     +(l+1)/(x^2)*hl(l,x)*(1-1/n^2)-(l+1)/x*(hl(l-1,x)-(l+1)/x*hl(l,x))*(1-1/n^2);

G = hl(l-2, x)+hl(l-1,x)*(-jl(l-1, n*x)/(n*jl(l,n*x))-l/x*(2-1/n^2))+...
    hl(l,x)*(-jl(l-2,n*x)/jl(l,n*x)+l*jl(l-1, n*x)/(n*x*jl(l, n*x))+(jl(l-1, n*x)/jl(l, n*x))^2+...
    (1-1/n^2)*l*(l+2)/x^2);

end

