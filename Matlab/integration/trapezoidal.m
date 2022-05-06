function val = trapezoidal(func,a,b,n)
% val = trapezoidal(FUNC,A,B,N)
%
% val is the approximation of \int_a^b func(x).dx using the trapezoidal rule
% func	-- function to integrate
% a	-- lower limit of integration interval
% b	-- upper limit of integration interval
% n	-- number of sub-intervals (= (number of evaluation points))
%
% func must be able to compute a vector of values given a vector of
% inputs (for efficiency)

val = func(a)+func(b);
h = (b-a)/n;
xs = (1:(n-1))*h + a;
val = val + 2*sum(func(xs));
val = val*h/2;
