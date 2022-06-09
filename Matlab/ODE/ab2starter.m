function [t1,y1] = ab2starter(f,y0,t0,h,m)
% function [t1,y1] = ab2starter(f,y0,t0,h,m)
%
% Given y(t0) = y0 and the step-size h > 0,
% we estimate y(t0+h) using a step-doubling method
% starting with an estimate for y(t0+h0) using
% Euler's method with h0 = 2^(-m)*h, and then
% using the 2nd order Adams-Bashforth method
% to double the step size.

h0 = h*2^(-m);
% Start with Euler's method with step-size h0
y1 = y0 + h0*f(t0,y0); % $y1 \approx y(t0+h0)$
% Now, use the step-size doubling technique m times,
% using the AB2 method:
for k = 1:m
    y2 = y1 + h0*((3/2)*f(t0+h0,y1) - (1/2)*f(t0,y0));
    % $ y2 \approx y(t0+2*h0) $
    y1 = y2;
    h0 = 2*h0;
    % $ y1 \approx y(t0+h0) $
end
t1 = t0 + h0;

    