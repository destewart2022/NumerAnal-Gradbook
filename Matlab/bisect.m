function [x,iter,xs] = bisect(func,a,b,eps,trace)
% function [x,iter,xs] = bisect(func,a,b,eps,trace)
%   Solves func(X) = 0 by bisection
%   Needs func(a) & func(b) to have opposite signs
%   eps is the tolerance for |b-a|
%   If trace is non-zero, trace execution
%   iter is the number of iterations used

if nargin < 5
    trace = 0;l
end
fa = func(a);
fb = func(b);
if ( fa == 0.0 )
  x = a;
  return
elseif ( fb == 0.0 )
  x = b;
  return
elseif ( (fa > 0.0 & fb > 0.0) | (fa < 0.0 & fb < 0.0 ) )
  fprintf('bisect: Error: signs of func(a) and func(b) are the same\n');
  return
end

xs = a;

iter = 0;
while abs(b-a) > eps
  if trace > 0
    fprintf('Iter # %d: a = %20.15g, f(a) = %10.6g, b = %20.15g, f(b) = %10.6g\n', ...
	iter,a,fa,b,fb);
  end
  c = a + 0.5*(b-a);
  fc = func(c);
  if ( fc == 0.0 )
    x = c;
    return
  end
  if ( (fa > 0.0 & fc > 0.0) | (fa < 0.0 & fc < 0.0) )
    % new interval is [c,b]
    a = c;
    fa = fc;
  else
    % new interval is [a,c]
    b = c;
    fb = fc;
  end
  iter = iter+1;
  if trace
    xs = [xs,c];
  end
end
x = a+0.5*(b-a);
