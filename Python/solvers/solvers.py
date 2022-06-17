# Various solvers in Python

import numpy as np

def bisect(f,a,b,eps,trace=0):
    """
    bisect(f,a,b,eps) applies bisection algorithm to approximately solve f(x) = 0
    The method terminates when |a-b| < eps and returns x = (a+b)/2.
    Prints trace of algorithm if trace > 0.
    Returns the pair (x,xs) 
    """
    fa = f(a);
    fb = f(b);
    if ( fa == 0.0 ):
        return a
    elif ( fb == 0.0 ):
        return b
    elif ( (fa > 0.0 and fb > 0.0) or (fa < 0.0 and fb < 0.0 ) ):
        print('bisect: Error: signs of f(a) and f(b) are the same\n')
        return null

    xs = [a]

    itern = 0
    while np.abs(b-a) > eps:
        if trace > 0:
            print(f"Iter # {itern}: a = {a}, f(a) = {fa}, b = {b}, f(b) = {fb}")
        c = a + 0.5*(b-a)
        fc = f(c)
        if ( fc == 0.0 ):
            return c
        if ( (fa > 0.0 and fc > 0.0) or (fa < 0.0 and fc < 0.0) ):
            # new interval is [c,b]
            a = c
            fa = fc
        else:
            # new interval is [a,c]
            b = c
            fb = fc
        itern += 1
        xs.append(c)
    x = a+0.5*(b-a);
    return (x,xs)

def secant(f,x0,x1,eps,maxit=100000,trace=0):
    """secant(f,x0,x1,eps,maxit=100000,trace=0) applies the secant method to solving f(x) = 0.
    Solves f(x) ~= 0 within accracy eps.
    Terminates when |f(x)| < eps or maxit iterations have been reached.
    """
    xs = [x0, x1]
    fx0 = f(x0)
    fx1 = f(x1)
    fvals = [fx0, fx1]
    itern = 0
    while abs(fx1) >= eps and itern < maxit:
        x_new = x1 - fx1/((fx1-fx0)/(x1-x0))
        if trace > 0:
            print(f"Iter {itern}: x0 = {x0}, x1 = {x1}, f(x0)={fx0}, f(x1)={fx1}")
        xs.append(x_new)
        fval_new = f(x_new)
        fvals.append(fval_new)
        # update x's
        x0 = x1
        x1 = x_new
        x = x1
        # and the function values
        fx0 = fx1
        fx1 = fval_new
        itern += 1
    return (x,xs,fvals)

def regfalsi(func,a,b,eps,maxit=100000,trace=0):
    """regfalsi(func,a,b,eps,maxit=100000,trace=0) uses the Regula Falsi method for solving f(x) = 0.
    Solves f(x) ~== 0 within accuracy eps.
    Terminates when |f(x)| < eps or maxit iterations have been reached.
    """
    fa = func(a)
    fb = func(b)
    if ( fa == 0.0 ):
        return a
    elif ( fb == 0.0 ):
        return b
    elif ( (fa > 0.0 and fb > 0.0) or (fa < 0.0 and fb < 0.0 ) ):
        print("regfalsi: Error: signs of func(a) and func(b) are the same");
        return

    xs = [a]

    iter = 0;
    while np.abs(b-a) > eps:
        if trace > 0:
            print(f"Iter # {iter}: a = {a}, f(a) = {fa}, b = {b}, f(b) = {fb}")
        c = a - (a-b)*fa/(fa-fb) # secant choice for next guess
        fc = func(c)
        if ( fc == 0.0 ):
            return c, xs
        if ( (fa > 0.0 and fc > 0.0) or (fa < 0.0 and fc < 0.0) ):
             # new interval is [c,b]
            a = c
            fa = fc
        else:
            # new interval is [a,c]
            b = c
            fb = fc
        iter = iter+1
        if trace:
            xs.append(c)
    x = c
    return (x, xs)

import numbers

def isnumber(x):
    """isnumber(x) returns a boolean according to whether x is a number (True) or is not (False)"""
    return isinstance(x,numbers.Number)

def newton(f,df,x,eps,maxit=10000,trace=0):
    """newton(f,df,x,eps,maxit=10000,trace=0) applies the (multivariate) Newton method to solving f(x) = 0.
    The method terminates when the norm of f(x) is less than eps or if the number of iterations reaches maxit.
    Returns the triple of (x,xs,fvals) where f(x) ~= 0, xs is the list of x values, and fvals[i] == f(xs[i]).
    """
    xs = [x];
    fval = f(x)
    fvals = [fval]
    itern = 0
    while np.linalg.norm(fval) >= eps and itern < maxit:
        dfval = df(x)
        if trace > 0:
            print(f"Iter {itern}: x = {x}, f(x) = {fval}, df(x) = {dfval}")
        if isnumber(dfval):
            x = x - fval / dfval
        else: # assume numpy matrix
            x = x - np.linalg.solve(dfval,fval)
        xs.append(x)
        fval = f(x)
        fvals.append(fval);
        itern += 1
    return (x,xs,fvals)

def gnewton(f,df,x,eps,maxit=10000,trace=0):
    """gnewton(f,df,x,eps,maxit=10000,trace=0) applies a guarded (multivariate) Newton method to solving f(x) = 0. The new x is given as x + s*d where x+d is the regular Newton update, where s > 0 is the largest negative power of two where ||f(x+s*d)|| <= (1-s/2)*||f(x)||.
    The method terminates when the norm of f(x) is less than eps or if the number of iterations reaches maxit.
    Returns the triple of (x,xs,fvals) where f(x) ~= 0, xs is the list of x values, and fvals[i] == f(xs[i]).
    """
    xs = [x];
    fval = f(x)
    fvals = [fval]
    itern = 0
    while np.linalg.norm(fval) >= eps and itern < maxit:
        dfval = df(x)
        if trace > 0:
            print(f"Iter {itern}: x = {x}, f(x) = {fval}, df(x) = {dfval}")
        if isnumber(dfval):
            d = - fval / dfval
        else: # assume numpy matrix
            d = - np.linalg.solve(dfval,fval)
        s = 1.0
        newfval = f(x+s*d)
        while np.linalg.norm(newfval) > (1-s/2)*np.linalg.norm(fval):
            s = s/2.0
            newfval = f(x+s*d)
            itern += 1
        x = x + s*d
        xs.append(x)
        fval = f(x)
        fvals.append(fval);
        if trace > 0:
            print(x,"  ",fval)
        itern += 1
    return (x,xs,fvals)

