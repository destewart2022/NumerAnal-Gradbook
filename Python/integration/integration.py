# Integration methods in Python

import numpy as np

def rectangle(f,a,b,n):
    """rectangle(f,a,b,n) estimates the integral of f(x) over [a,b] using n sub-intervals and the rectangle rule.
    
    Note that this uses exactly n function evaluations of f.
    """
    h = (b-a)/n
    return h*sum(f(a+i*h) for i in range(n))

def trapezoidal(f,a,b,n):
    """trapezoidal(f,a,b,n) estimates the integral of f(x) over [a,b] using n sub-intervals and the trapezoidal rule.
    
    Note that this uses n+1 function evaluations of f."""
    h = (b-a)/n
    return h*((f(a)+f(b))/2+sum(f(a+i*h) for i in range(1,n)))

def simpson(f,a,b,n):
    """simpson(f,a,b,n) estimates the integral of f(x) over [a,b] using n sub-intervals and Simpson's rule.
    
    Note that this uses 2n+1 function evaluations of f."""
    h = (b-a)/n
    return (h/6)*(f(a)+f(b) +
                  4*sum(f(a+(i+1/2)*h) for i in range(0,n)) + 
                  2*sum(f(a+      i*h) for i in range(1,n)))

def midpt(f,a,b,n):
    """midpt(f,a,b,n) estimates the integral of f(x) over [a,b] using n sub-intervals and the midpoint rule.
    
    Note that this uses n function evaluations of f."""
    h = (b-a)/n
    return h*sum(f(a+(i+1/2)*h) for i in range(n))

def gen_rule(f,ws,xs,a,b,n):
    """gen_rule(f,a,b,n) estimates the integral of f(x) over [a,b] using n sub-intervals and on each sub-interval uses sum of ws[i]*f(a+(j+xs[i])*h) over i where h is the width of each sub-interval; j is the sub-interval number.
    
    The pair (ws,xs) should define an integration method over the interval [0,1].
    
    Note that this uses n*len(ws) function evaluations of f."""
    h = (b-a)/n
    return h*sum(sum(ws[i]*f(a+h*(j+xs[i])) for i in range(len(ws))) for j in range(n))