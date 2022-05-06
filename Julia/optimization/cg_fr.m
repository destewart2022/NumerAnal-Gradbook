function [x,fval,g,nfe,nge,xs] = ...
    cg_fr(func,x,tol,c1,c2,it_limit,trace)
% [X,FVAL,GVAL,NFE,NGE,XS] = ...
%	cg_fr(FUNC,X,TOL,C1,C2,IT_LIMIT,TRACE)
% Fletcher-Reeves conjugate gradient method (without restart)
% Uses Wolfe line-search conditions (parameters c1 and c2)
%
% FUNC is the name of the function to minimize: FUNC(X)
% DFUNC is the name of the gradient function: DFUNC(X)
% X is the current position
% TOL is the stopping tolerance: ||gradient|| < tol
% C1 and C2 are the Wolfe line-search parameters
%	(See wolfe_search2.m.)
% IT_LIMIT is the iteration limit
% If TRACE is non-zero, then print out information about process
%
% Returns
%	X	-- final "minimizer"
%	FVAL	-- FUNC(X), function value at X
%	GVAL	-- gradient value at X
%	NFE	-- # function eval'ns
%	NGE	-- # gradient eval'ns
%	XS	-- list of iterates (only used if trace ~= 0)

% Initialization
if nargin <= 6
  trace = 0
end

xs = [];

[fval,g] = func(x);
nfe = 1;
nge = 1;
p = -g;
g_sqr = g'*g;
for k = 0:it_limit
  if trace > 0
    xs = [xs x];
  end
  if trace > 1
    fprintf('cg_fr: cos(theta) = %g\n', -p'*g/(norm(p,2)*norm(g,2)))
    fprintf('cg_fr: ||p||/||g|| = %g\n', norm(p,2)/norm(g,2))
    pause
  end
  [alpha,fval,new_g,nfe_ws,nge_ws] = ...
      wolfe_search2(func,x,p,1.0,c1,c2,fval,g,max(trace-1,0));
  nfe = nfe + nfe_ws;
  nge = nge + nge_ws;
  if trace ~= 0
    fprintf('cg_fr: iter # %d, alpha = %g, function value = %g\n', ...
	k, alpha, fval);
  end
  x = x + alpha*p;
  new_g_sqr = new_g'*new_g;
  beta = new_g_sqr/g_sqr;
  g_sqr = new_g_sqr;
  g = new_g;
  if norm(g,2) < tol
    return
  end
  % Restart code (if desired)
  % if k % n == n-1
  %   beta = 0;
  % end
  p = -g + beta*p;
end
fprintf('Have used %d iterations; returning without finding solution\n', ...
    it_limit);

