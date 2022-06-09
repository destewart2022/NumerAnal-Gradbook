function [y,t] = ab2(func,y0,y1,t0,h,n)
% [Y,T] = ab2(FUNC,Y0,Y1,T0,H,N)
% Compute solution of ODE y' = func(t,y), y(t0) = y0 using the 2nd order
% Adams-Bashforth method:
%   y_{n+1} = y_n + h*(1.5*y_n' - 0.5*y_{n-1}')
% Y1 is approximation to y(t0+h)
% Computes solution with step-size H for N steps.
% Solution stored in y so that y(t(i)) ~= y(:,i)

y  = zeros(length(y0),n+1);
fy = zeros(length(y0),n+1);
t  = t0+h*(0:n);
y(:,1) = y0;
y(:,2) = y1;

fy(:,1) = feval(func,t0,y0);
for i = 2:n
  fy(:,i) = feval(func,t(i),y(:,i));
  y(:,i+1) = y(:,i) + h*(1.5*fy(:,i)-0.5*fy(:,i-1));
end
