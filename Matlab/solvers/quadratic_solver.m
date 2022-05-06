function [xp,xm] = quadratic_solver(a,b,c)
% function [xp,xm] = quadratic_solver(a,b,c)
%
% Returns both solutions to a*x^2+b*x+c = 0
if b < 0
    xp = (-b+sqrt(b^2-4*a*c))/(2*a);
    xm = (c/a)/xp;
else
    xm = (-b-sqrt(b^2-4*a*c))/(2*a);
    xp = (c/a)/xm;
end

end

