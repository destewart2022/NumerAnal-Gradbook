function rtable = romberg(f,a,b,n)
% function rtable = romberg(f,a,b,n)
%
% Compute the Romberg table for integrating f from a to b
% using an n x n Romberg table.
% This will require evaluating f at 2^n points.

rtable = zeros(n,n);

% Evaluate trapezoidal integrals

% Evaluate simple trapezoidal rule
rtable(1,1) = 0.5*(b-a)*(feval(f,a)+feval(f,b));

% Now evaluate at the other points
for k = 2:n
  pow2km1 = 2^(k-1);
  h = (b-a)/pow2km1;
  rtable(k,1) = 0.5*rtable(k-1,1) + h*sum(feval(f,a+h*(1:2:pow2km1)));
end

% Now construct the Romberg table
for k = 2:n
  for i = k:n
    rtable(i,k) = (4^(k-1)*rtable(i,k-1)-rtable(i-1,k-1))/(4^(k-1)-1);
  end
end
