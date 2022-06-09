function [y,fy] = ode_solver(f,t,w,y,alpha,eps,params)
% function [y,fy] = ode_solver(f,t,w,y,alpha,eps,params)
%
% Interface for solvers to pass to implicit ODE routines
% for solving
%    y = w + alpha*f(t,y)
% so that the returned y satisfies
%    ||y-(w + alpha*f(t,y))|| <= eps.
% On return, fy = f(t,y).
% This function is empty as it merely describes the 
% required interface.

this will force a stop;  % force error
