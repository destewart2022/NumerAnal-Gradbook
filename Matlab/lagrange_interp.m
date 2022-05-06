function Lvals = lagrange_interp(i,x,xs)
% function Lvals = lagrange_interp(i,x,xs)
%
% Implements Lagrange interpolation polynomial L_i(x)
% at x (or its entries) with interpolation points xs(i)
Lvals = ones(size(x));
for j = [1:i-1,i+1:length(xs)]
    Lvals = Lvals .* ((x-xs(j))/(xs(i)-xs(j)));
end
