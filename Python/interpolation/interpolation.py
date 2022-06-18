# Interpolation and approximation codes in Python

import numpy as np

def divdif(xs,ys):
    """divdif(xs,ys) computes divided difference table:
    
    returns [f(x0), f[x0,x1], f[x0,x1,x2], ...] where f[x0,x1,...,xn] is the divided 
    difference of f with f(xj) = yj."""
    n = len(xs)
    dd = np.copy(ys)
    for k in range(n):
        # compute (k+1)'th order divided differences
        # print(f"k = {k}")
        for j in range(n-1,k,-1):
            # print(f"j = {j}")
            dd[j] = (dd[j]-dd[j-1]) / (xs[j] - xs[j-k-1])
    return dd

def dd_interp(xs,dd,t):
    """dd_interp(xs,dd,t) return the value of the 
    polynomial interpolant given by interpolation points xs and divided difference table dd at t.
    """
    n = len(xs) - 1
    pval = dd[n]
    for i in range(n-1,-1,-1):
        pval = dd[i] + (t-xs[i])*pval
    return pval

