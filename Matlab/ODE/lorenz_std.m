function rhs = lorenz_std(t,u,params)
% function rhs = lorenz_std(t,u,params)
%
% Computes right-hand side of Lorenz' famous
% "chaotic" ODE:
%   dx/dt = -s*x+s*y
%   dy/dt = -x*z+r*x-y
%   dz/dt = x*y-b*z
%
% u = [x; y; z]; params contains r, s, and b
%
% Standard values: s = 10; r = 28; b = 8/3
  rhs = zeros(3,1);
  rhs(1) = params.s*(u(2)-u(1));
  rhs(2) = u(1)*(params.r-u(3))-u(2);
  rhs(3) = u(1)*u(2)-params.b*u(3);
  