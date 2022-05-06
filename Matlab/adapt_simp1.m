function [val,n_f_evals] = adapt_simp1(func,a,b,f0,f12,f1,tol,trace)
% function val = adapt_simp1(func,a,b,f0,f12,f1,tol,trace)
%
% Recursive routine for adaptive Simpson's rule.
% The error is estimated by computing the integral
% over [a,b] using Simpson's rule, and composite Simpson's
% rule using two subintervals.  The error estimate
% is then (1/15)*(the difference in the integrals).
% If the error estimate is more than tol*(b-a), then
% adapt_simp1 is called recursively on the sub-divided intervals
% [a, (a+b)/2] and [(a+b)/2, b].
%
% Inputs:
%    func	= name of function to integrate
%    a, b	= endpoints of integration interval
%    f0, f12, f1 = values func(a), func((a+b)/2), and f(b) respectively
%    tol	= error tolerance divided by the interval length.
%
% Note: f0, f12, f1 are passed to avoid re-evaluating the function
% func unnecessarily.

% Compute integral estimates via Simpson's rule:
int1 = (b-a)/6*(f0+4*f12+f1);
f14 = func(a+.25*(b-a));
f34 = func(a+.75*(b-a));
int2 = (b-a)/12*(f0+4*f14+2*f12+4*f34+f1);
n_f_evals = 2;

% Estimate error
err_est = (int2-int1)/15;
if trace > 0
    fprintf('adapt_simp1: error estimate for [%g,%g] = %g\n', ...
        a,b,err_est)
end

if abs(err_est) > tol
  % Recursively call adapt_simp1 on sub-intervals
  [val1,n_f_evals1] = adapt_simp1(func,a,(a+b)/2,f0,f14,f12,tol);
  [val2,n_f_evals2] = adapt_simp1(func,(a+b)/2,b,f12,f34,f1,tol);
  val = val1 + val2;
  n_f_evals = n_f_evals + n_f_evals1 + n_f_evals2;
else
  % Return the most accurate estimate of the integral
  val = int2;
end
