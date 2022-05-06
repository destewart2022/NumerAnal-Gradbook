function pvals = dd_interp(xs,dd,t)
% PVALS = dd_interp(XS,DD,T)
% Evaluates interpolations polynomial defined
%	by interpolation points XS and
%	the divided difference table DD
% Evaluate the interpolating polynomial at T
% If T is a vector or matrix, then the interpolating
% polynomial is applied componentwise.
%
% See also divdif.m.
%
  n = length(xs) - 1;
  pvals = ones(size(t))*dd(n+1);
  for i = n-1:-1:0
    pvals = dd(i+1) + (t-xs(i+1)).*pvals;
  end
  
