function [y,t] = bdf3(func,solver,y0,y1,y2,t0,h,n,eps,params)
% [y,t] = bdf3(func,solver,y0,y1,y2,t0,h,n,eps,params)
% Compute solution of ODE y' = func(t,y), y(t0) = y0 using 
% the 3rd order BDF method:
%   y_{n+1} = (18/11)*y_n - (9/11)*y_{n-1} + (2/11)*y_{n-2}
%                             + (h*6/11)*f(t_{n+1},y_{n+1})
% Computes solution with step-size h for n steps.
% Solution stored in y so that y(t(i)) ~= y(:,i)
%
% Implicit equations solved by solver(...), which has the same
% interface as ode_solver() (see ode_solver.m).
%
% y0 is initial value.
% y1 is estimate for y(t1) = y(t0+h)
% y2 is estimate for y(t2) = y(t0+2*h)

y  = zeros(length(y0),n+1);
t  = t0 + h*(0:n);
y(:,1) = y0;
y(:,2) = y1;
y(:,3) = y2;

t_val =  t0+2*h;
y_old2 = y0;
y_old1 = y1;
y_val  = y2;

fy_old1 = feval(func,t0+h,  y_old1);
fy_val  = feval(func,t_val, y_val);
for i = 3:n
  % Predictor: AB-2
  y_new = y_val + h*(1.5*fy_val - 0.5*fy_old1);
  % Now use solver to solve implicit equations
  w = (18/11)*y_val - (9/11)*y_old1 + (2/11)*y_old2;
  [y_new,fy_new] = feval(solver,func,t_val+h,w,y_new,h*6/11,eps,params);
  norm_solver_error = norm(y_new-((18/11)*y_val - (9/11)*y_old1 + (2/11)*y_old2 ...
      + h*6/11*fy_new))
  % Update fy_... and y_... values
  fy_old2 = fy_old1;
  fy_old1 = fy_val;
  fy_val  = fy_new;
  y_old2  =  y_old1;
  y_old1  =  y_val;
  y_val   =  y_new;
  
  % Store new value
  y(:,i+1) = y_val;
  t_val = t_val + h;
end
