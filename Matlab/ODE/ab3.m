function [y,t] = ab3(func,y0,y1,y2,t0,h,n)
% [Y,T] = ab3(func,y0,y1,y2,t0,h,n)
% Compute solution of ODE y' = func(t,y), y(t0) = y0 using the 2nd order
% Adams-Bashforth method:
%   y_{n+1} = y_n + h/12*(23*y_n' - 16*y_{n-1}' + 5*y_{n-2}')
% Y1 is approximation to y(t0+h)
% Y2 is approximation to y(t0+2*h)
% Computes solution with step-size h for n steps.
% Solution stored in y so that y(t(i)) ~= y(:,i)

y  = zeros(length(y0),n+1);
fy = zeros(length(y0),n+1);
t  = t0+h*(0:n);
y(:,1) = y0;
y(:,2) = y1;
y(:,3) = y2;

fy(:,1) = feval(func,t0,y0);
fy(:,2) = feval(func,t0+h,y1);
for i = 3:n
  fy(:,i) = feval(func,t(i),y(i));
  y(:,i+1) = y(:,i) + (h/12)*(23*fy(:,i)-16*fy(:,i-1)+5*fy(:,i-2));
end
