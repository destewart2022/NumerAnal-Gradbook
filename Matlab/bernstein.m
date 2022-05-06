function yvals = bernstein(xvals,k,n)
% yvals = bernstein(xvals,k,n)
%
% yvals(i) = phi_{k,n}(xvals(i))
% where phi_{k,n}(x) = (n choose k)x^k(1-x)^(n-k)
%
% This routine uses logarithms to compute this
% for 0 < xvals < 1
ln_yvals = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+ ...
    k*log(xvals)+(n-k)*log(1-xvals);
yvals = exp(ln_yvals);
