function [x,xs,fvals] = gnewton(func,dfunc,x,eps,trace)
% [x,xs,fvals] = gnewton(func,dfunc,x,eps,trace)
%
% Uses a guarded form of Newton's method to solve
%	func(x) = 0
% Approximate solution returned as x.
% A history of the iterates is returned as xs and the function
% values as fvals.
% The input x is the initial guess.
% The function is given as a Matlab function with the name in func
% and its derivative is another Matlab function with the name in dfunc.
%
% We compute dx = -fval/dfval = -f(x)/f'(x).
% We then find and alpha > 0 where |f(x+alpha.dx)| < |f(x)|
% starting at alpha = 1 and halving until we find a suitable value of alpha.
% The new value of x is x+alpha*dx.
if nargin < 5
    trace = 0;
end
xs = [x];
fval = feval(func,x);
fvals = [fval];
while norm(fval) >= eps 
  dfval = dfunc(x);
  dx = - dfval\fval;
  alpha = 1;
  x_new = x + alpha*dx;
  fval_new = func(x_new);
  while norm(fval_new) >= (1-alpha/2)*norm(fval) & alpha >= 1e-12
    alpha = 0.5*alpha;
    x_new = x + alpha*dx;
    fval_new = func(x_new);
  end
  if trace > 0
      x_new
      fval_new
      alpha
  end
  if alpha < 1e-12
      break
  end
  fval = fval_new;
  fvals = [fvals fval];
  x = x_new;
  xs = [xs x];
end
