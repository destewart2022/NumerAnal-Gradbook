function x = newton(func,dfunc,x,eps,trace)
% newton(func,dfunc,x,eps,trace)
% Performs one step of Newton's method from X.
%   func   -- function to solve func(x) = 0
%   dfunc  -- derivative/Jacobian function 
%   eps    -- stopping tolerance
%   trace  -- set to non-zero for trace of execution
% Returns the new value x
  fx = func(x);
  iter = 0;
  while norm(fx,2) >= eps
    dfx = dfunc(x);
    if trace
      fprintf('Iter # %d: [x,fx] = ',iter);
      [x,fx]
      norm_fx = norm(fx)
    end
    x = x - dfx\fx;
    fx = func(x);
    iter = iter + 1;
  end
  if trace
    fprintf('Iter # %d: [x,fx] = ',iter);
    [x,fx]
    norm_fx = norm(fx)
  end
end % function

