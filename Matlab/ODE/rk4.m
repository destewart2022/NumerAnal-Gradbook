function [tlist,ylist] = rk4(f,y0,t0,h,n)
% function [tlist,ylist] = rk4(f,y0,t0,h,n)
% Carry out n steps of the standard 4th order 
% RK method for dy/dt = f(t,y)
% with step size h and y(t0) = y0
ylist = zeros(length(y0),n+1);
tlist = zeros(1,n+1);
tlist(1) = t0;
ylist(:,1) = y0;
y = y0;
for i = 1:n
  v1 = feval(f,t0+(i-1)*h,y);
  v2 = feval(f,t0+(i-.5)*h,y+(h/2)*v1);
  v3 = feval(f,t0+(i-.5)*h,y+(h/2)*v2);
  v4 = feval(f,t0+i*h,y+h*v3);
  y = y + (h/6)*(v1+2*v2+2*v3+v4);
  ylist(:,i+1) = y;
  tlist(i+1) = t0+i*h;
end

