function integral = trap(a,b,n,f)
% function integral = trap(a,b,n,f)
%
% This uses the trapezoidal rule to compute
% the integral of f(x).dx from x = a to x = b.
% This uses n pieces.

h = (b-a)/n;
s = 0.5*(f(a)+f(b));
s = s + f(a+(1:n-1)*h);
integral = s*h;
