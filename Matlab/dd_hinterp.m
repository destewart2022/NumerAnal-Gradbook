function pvals = dd_hinterp(xs,dd,t)
% pvals = dd_hinterp(xs,dd,t)
% Evaluates interpolations polynomial defined
%	by interpolation points XS and
%	the divided difference table DD
% Evaluate the interpolating polynomial at T
% If T is a vector or matrix, then the interpolating
% polynomial is applied componentwise.
%
% See also dd_interp.m and dd_hermite.m.
%
  new_xs = zeros(1,2*length(xs));
  new_xs(2*(1:length(xs)))   = xs;
  new_xs(2*(1:length(xs))-1) = xs;
  pvals = dd_interp(new_xs,dd,t);

  