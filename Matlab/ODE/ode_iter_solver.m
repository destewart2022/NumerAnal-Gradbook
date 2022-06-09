function [y,fy] = ode_iter_solver(f,t,w,y,alpha,eps,params)
% function y = ode_iter_solver(f,t,w,y,alpha,eps,params)
% 
% Simple iterative solver for ode_solver() interface.
fy = feval(f,t,y,params);
while norm(y - (w+alpha*fy)) >= eps
    y = w + alpha*fy;
    fy = feval(f,t,y,params);
end
