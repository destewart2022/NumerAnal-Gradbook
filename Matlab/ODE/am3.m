function [y,t] = am3(func,solver,y0,y1,t0,h,n,eps,params)
% [y,t] = am3(func,solver,y0,y1,t0,h,n,eps,params)
% Compute solution of ODE y' = func(t,y,params), y(t0) = y0 using the 3rd order
% Adams-Moulton method:
%   y_{n+1} - y_n = (h/12)*(5*y_{n+1}' + 8*y_n' - y_{n-1}')
% Computes solution with step-size h for n steps.
% Solution stored in y so that Y(t(i)) ~= y(:,i)
%
% Implicit equations solved by solver(...), which has the same
% interface as ode_solver() (see ode_solver.m).
%
% y0 is initial value.
% y1 is estimate for y(t1) = y(t0+h)
%
% Uses simple fixed-point iteration to solve implicit equations,
% until residual is no more than eps in 2-norm.
%
% Uses the 2nd order Adams--Bashforth method for the predictor.

y  = zeros(length(y0),n+1);
t  = t0 + h*(0:n);
y(:,1) = y0;
y(:,2) = y1;

fy_old2 = feval(func,t0  ,y0,params);
fy_old1 = feval(func,t0+h,y1,params);
t_val =  t0+h;
y_old1 = y0;
y_val  = y1;
for i = 2:n
  y_old2 = y_old1;
  y_old1 = y_val;
  % AB2 predictor
  y_val = y_old1 + h*(1.5*fy_old1 - 0.5*fy_old2);
  % Now use solver to solve implicit equations
  w = y_old1 + h*(2/3*fy_old1-1/12*fy_old2);
  [y_val,fy_val] = feval(solver,func,t_val+h,w,y_val,h*5/12,eps,params);
  %norm_solver_error = norm(y_val-(y_old1+h*(5/12*fy_val+2/3*fy_old1-1/12*fy_old2)))
  % Save fy_val for next iteration
  fy_old2 = fy_old1;
  fy_old1 = fy_val;
  % Store new value
  y(:,i+1) = y_val;
  t_val = t_val + h;
end
