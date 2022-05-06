function dd = divdif(xs,fs)
% DD = divdif(XS,FS)
% Construct divided difference table
%
%   XS = array of interpolation points
%   FS = array of (function) values to interpolate
% Returns DD, the divided difference table
%
% See also dd_interp.m.
%
  n = length(xs) - 1;
  dd = zeros(1,n+1);
  for i = 1:n+1
    dd(i) = fs(i);
  end
  
  for k = 1:n  
    % compute k'th order divided differences
    for j = n:-1:k
      dd(j+1) = (dd(j+1) - dd(j))/(xs(j+1) - xs(j-k+1));
    end
  end
  