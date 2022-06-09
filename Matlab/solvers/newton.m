function x = newton(func,dfunc,x,eps,trace)
% newton(FUNC,DFUNC,X,EPS,TRACE)
% Performs one step of Newton's method from X.
%   FUNC   -- function name
%   DFUNC  -- gradient function name
%   EPS    -- stopping tolerance
%   TRACE  -- set to non-zero for trace of execution
% Returns the new value X
  fx = feval(func,x);
  iter = 0;
  while norm(fx,2) >= eps
    dfx = feval(dfunc,x);
    if trace
      fprintf('Iter # %d: [x,fx] = ',iter);
      [x,fx]
      norm_fx = norm(fx)
    end
    x = x - dfx\fx;
    fx = feval(func,x);
    iter = iter + 1;
  end
  if trace
    fprintf('Iter # %d: [x,fx] = ',iter);
    [x,fx]
    norm_fx = norm(fx)
  end
end % function

