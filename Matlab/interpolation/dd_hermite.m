function dd = dd_hermite(xs,fs,dfs)
% dd = dd_hermite(xs,fs,dfs)
% Construct divided difference table for Hermite interpolant
%
%   xs = array of interpolation points
%   fs = array of (function) values to interpolate
%   dfs= array of derivative values to interpolate
% Returns dd, the divided difference table
%
% See also dd_interp.m.
%
  n = length(xs) - 1;
  dd = zeros(1,2*n+2);
  dd(1) = fs(1);
  dd(2) = dfs(1);
  for i = 1:n
    dd(2*i+1) = (fs(i+1)-fs(i))/(xs(i+1)-xs(i));
    dd(2*i+2) = dfs(i+1);
  end
  
  new_xs = zeros(1,2*length(xs));
  new_xs(2*(1:length(xs)))   = xs;
  new_xs(2*(1:length(xs))-1) = xs;
  for i = 2:(2*n+1)
    for j = (2*n+1):-1:i
      dd(j+1) = (dd(j+1) - dd(j)) / (new_xs(j+1) - new_xs(j-i+1));
    end
  end
  
