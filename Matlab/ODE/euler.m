function [t,y] = euler(func,y0,t0,h,n)
% function [t,y] = euler(func,y0,t0,h,n)
%
% Compute solution of ODE y' = func(t,y), y(t0) = y0 using Euler's
% method:
%   $y_{n+1} = y_n + h*func(t_n,y_n)$
% Computes solution with step-size h for h steps.
% Solution stored in y so that y_i = y(:,i); $y_i \approx y(t(i))$,
% $t(i)=t_0+i\,h$.
y = zeros(length(y0),n+1);
t = t0 + h*(0:n);
y(:,1) = y0;

for i = 1:n
  fy = feval(func,t(i),y(:,i)); % or fy = func(t(i),y(:,i))
  y(:,i+1) = y(:,i) + h*fy;
end
