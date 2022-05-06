function pval = horner(coeffs,x)
% function pval = horner(coefs,x)
%
% Evaluates the polynomial sum(i=0..n) coefs(i+1)*x^i efficiently
% via Horner's nested method
n = length(coeffs)-1;
pval = coeffs(n+1)*ones(size(x));
for k = n-1:-1:0
	pval = coeffs(k+1) + x.*pval;
end

