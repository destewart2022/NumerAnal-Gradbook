function [val,n_f_evals] = adapt_simp(func,a,b,tol,trace)
% function [val,n_f_evals] = adapt_simp(func,a,b,tol,trace)
%
% Adaptive Simpson's rule to compute 
% 
%     integral_a^b func(x).dx
%
% to within an error of tol*(b-a).
% Calls recursive routine adapt_simp1.
%
% Inputs:
%    func	= name of function to integrate
%    a, b	= endpoints of integration interval
%    tol	= error tolerance divided by the interval length.
%
% Outputs:
%    val	= computed value of the integral
%    n_f_evals	= number of function evaluations used
if nargin < 5
    trace = 0;
end
f0 = func(a);
f1 = func(b);
f12 = func((a+b)/2);
n_f_evals = 3;
[val,n_f_evals1] = adapt_simp1(func,a,b,f0,f12,f1,tol,trace);
n_f_evals = n_f_evals + n_f_evals1;
