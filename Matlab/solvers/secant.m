function [x,xs,fvals] = secant(func,x1,x2,eps,trace)
% [x,xs,fvals] = secant(func,x1,x2,eps,trace)
%
% Uses the secant method to solve
%	func(x) = 0
% Approximate solution returned as x.
% A history of the iterates is returned as xs and the function
% values as fvals.
% The inputs x1 and x2 are the initial guess.
% The function is given as a Matlab function with the name in func
% and its derivative is another Matlab function with the name in dfunc.

if nargin < 5
    trace = 0;
end
x = x2;
xs = [x1, x2];
fval1 = func(x1);
fval2 = func(x2);
fvals = [fval1 fval2];
while abs(fval2) >= eps && fval1 ~= fval2
  x_new = x2 - fval2/((fval2-fval1)/(x2-x1));
  xs = [xs x_new];
  fval_new = feval(func,x_new);
  if trace > 0
      x_new
      fval_new
  end
  fvals = [fvals fval2];
  % update x's
  x1 = x2;
  x2 = x_new;
  x = x2;
  % and the function values
  fval1 = fval2;
  fval2 = fval_new;
end
