function val = simpson(func,a,b,n)
% VAL = simpson(FUNC,A,B,N)
%
% VAL is the approximation of \int_a^b func(x).dx using Simpson's rule
% FUNC	-- function to integrate
% A	-- lower limit of integration interval
% B	-- upper limit of integration interval
% N	-- number of sub-intervals (= (number of evaluation points)/2)
%
% FUNC must be able to compute a vector of values given a vector of
% inputs (for efficiency)


%  val = (b-a)/6*(feval(func,a)+4*feval(func,(a+b)/2)+feval(func,b));

val = feval(func,a)+feval(func,b);
h = (b-a)/n;
xs = ((0:(n-1))+0.5)*h + a;
val = val + 4*sum(feval(func,xs));
xs = (1:(n-1))*h + a;
val = val + 2*sum(feval(func,xs));
val = val*h/6;
