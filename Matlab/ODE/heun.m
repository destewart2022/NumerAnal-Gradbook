function [t,y] = heun(func,y0,t0,h,n)
% [T,Y] = heun(FUNC,Y0,T0,H,N)
% Compute solution of ODE y' = func(t,y), y(T0) = Y0 using Heun's
% method:
%   z_{n+1} - y_n = H*f(t_n,y_n)
%   y_{n+1} - y_n = (H/2)*(f(t_n,y_n)+f(t_{n+1},z_{n+1})
% Computes solution with step-size H for N steps.
% Solution stored in Y so that y(t(i)) ~= Y(i)
y = zeros(length(y0),n+1);
t = t0 + h*(0:n);
y(:,1) = y0;

for i = 1:n
  fy = feval(func,t(i),y(:,i));
  z = y(:,i) + h*fy;
  fz = feval(func,t(i)+h,z);
  y(:,i+1) = y(:,i) + (h/2)*(fy+fz);
end
