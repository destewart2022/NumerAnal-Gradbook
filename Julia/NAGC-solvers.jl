"""
    NAGC Solvers module

Includes `bisect`, `fixedpt`, `secant`, `regfalsi`, `newton` and `gnewton` functions for solving ``f(x) = 0`` or ``g(x) = x`` .
"""

module NAGCSolvers

using LinearAlgebra

export bisect, fixedpt, regfalsi, secant, newton, gnewton



""" 
x,nfe = bisect(f::Function,a::T,b::T,eps::T;
    max_iter::Int=100, trace::Int=0, xlist::Union{Vector{T},Nothing}=nothing)

    Uses bisection method on interval ``[a,b]`` (or ``[b,a]`` if ``a > b``) to
    compute approximate solution `x` to `f(x) == 0`.

### Arguments
- `f(x)` function 
- `a` 1st endpoint of interval
- `b` 2nd endpoint of interval
- `eps` target accuracy
### Optional named arguments
- `max_iter` maximum number of iterations (default is 100)
- `trace` trace execution of algorithm (default is 0 for no tracing)
"""
function bisect(f::Function,a::T,b::T,eps::T;
        max_iter::Int=100, trace::Int=0, xlist::Union{Vector{T},Nothing}=nothing) where {T}
fa = f(a); fb = f(b);
nfe = 2

if sign(fa)*sign(fb) > 0
    println("f(a0) and f(b0) are of the same sign.  Stop!")
    return
end

c = (a+b)/2;
it_count = 0;
if trace > 0
    println(" iter\t  a\t  b\t  c\t  f(c)\t b-c")
end
while b-c > eps && it_count < max_iter
    it_count = it_count + 1
    fc = f(c)
    nfe = nfe + 1
    if trace > 0
        println(it_count,"\t", a,"\t", b,"\t", c,"\t", fc,"\t\t", b-c)
    end
    if ( sign(fb)*sign(fc) <= 0 )
        a = c
        fa = fc
    else
        b = c
        fb = fc
    end
    c = (a+b)/2
    if xlist != nothing; push!(xlist,c) end
end
return (c,nfe)
end 



"""
    fixedpt(g::Function,x::V,eps::T;
        max_iter::Int=100, trace::Int=0, xlist::Union{Vector{T},Nothing}=nothing) where {T}

Solves ``g(x)=x`` via the iteration ``x_{n+1}\\gets g(x_n)``.

Returns iterate `x` and number of function evaluations `nfe`.

### Arguments
-  `f(x)` is the function value at `x`,
-  `x` (input) is the initial guess,
-  `eps > 0` is the error tolerance: the method stops when `||f(x)|| < eps`.

### Optional named arguments
-  `xlist` is the list of the iterates (if `xlist != nothing`)
-  `trace=1` to trace the execuation of the algorithm
-  `max_iter=m` to limit the algorithm to no more than `m` iterations (default is 10000)

"""
function fixedpt(g::Function,x::V,eps::T;
        max_iter::Int=100, trace::Int=0, xlist::Union{Vector{V},Nothing}=nothing) where {T,V}
    gx = g(x)
    nfe = 1
    while norm(gx-x) >= eps && nfe <= max_iter
        x = gx
        gx = g(x)
        nfe = nfe + 1
        if trace > 0
            println("fixedpt: x = ",x)
        end
        if xlist != nothing
            push!(xlist,x) 
        end
    end
    return x, nfe
end # function


""" 
    x,nfe = secant((func::Function,x1::T,x2::T,eps::T;
        max_iter::Int=100, trace::Int=0, xlist::Union{Vector{T},Nothing}=nothing) where {T}

 Uses the secant method to solve
    `func(x) = 0`
 Approximate solution returned as `x`; `nfe` is the number of function evaluations.
### Arguments
- `func` function to solve for
- `x1` and `x2` are the initial guesses.
### Optional names arguments
- `xlist` a history of the iterates if `xlist != nothing`
- `max_iter` the maximum number of iterations (default is 100)
- `trace` for tracing the operation of the algorithm (default is zero for no tracing)
 
"""
function secant(func::Function,x1::T,x2::T,eps::T;
        max_iter::Int=100, trace::Int=0, xlist::Union{Vector{T},Nothing}=nothing) where {T}

x = x2
fval1 = func(x1)
fval2 = func(x2)
nfe = 2

while abs(fval2) >= eps
  if ! iszero(fval2 - fval1)
    x_new = x2 - fval2*(x2-x1)/(fval2-fval1)
  else
    return x2,nfe
  end
  if xlist != nothing; push!(xlist,x_new) end
  if trace > 0
    println("secant: new x = ",x,", after ",nfe," function evaluations")
  end
  fval_new = func(x_new)
  nfe = nfe + 1
  # update x's
  x1 = x2
  x2 = x_new
  x = x2
  # and the function values
  fval1 = fval2
  fval2 = fval_new
  if xlist != nothing; push!(xlist,x2) end
end
return x2,nfe
end # function


"""
    function newton(f::Function, df::Function, x::V, eps::Number; 
                  trace=0::Int, itmax=10000::Int, xlist::Union{Vector{V},Nothing}=nothing)

  Implements Newton's method for solving `f(x) == 0`.
### Arguments
-  `f(x)` is the function value at `x`,
-  `df(x)` is the Jacobian matrix of `f` at `x`,
-  `x` (input) is the initial guess,
-  `eps > 0` is the error tolerance: the method stops when `||f(x)|| < eps`.

### Optional named arguments
-  `xlist` is the list of the iterates (if `xlist != nothing`)
-  `trace=1` to trace the execuation of the algorithm
-  `itmax=m` to limit the algorithm to no more than `m` iterations (default is 10000)

  Returns solution, along with number of function and derivative (i.e., Jacobian) evaluations (`nfe`, `nde` resp.)

  Example: To solve ``x-\\cos(x) = 0`` with a tolerance of ``10^{-10}`` and starting at ``x=0``:
    ```x,nfe,nde = newton(x->x-cos(x),x->1+sin(x),0,1e-10)```
"""
function newton(f::Function, df::Function, x::V, eps::Number; 
        trace::Int=0, itmax::Int=10000, xlist::Union{Vector{V},Nothing} = nothing) where {V}
fx = f(x)
nfe = 1
nde = 0
iter = 0
while norm(fx) >= eps && iter <= itmax
  dfx = df(x)
  nde = nde + 1
  x = x - (dfx \ fx)
  fx = f(x)
  nfe = nfe + 1
  if trace > 1
    println("newton: iter = ", iter, ", x = ", x,", fx = ",fx,", dfx (old) = ",dfx)
  elseif trace > 0
    println("newton: iter = ", iter, ", ||f(x)|| = ",norm(fx))
  end
  if xlist != nothing
    push!(xlist,x)
  end
  iter = iter + 1
end
x,nfe,nde # returns solution, along with number of function and derivative (i.e., Jacobian) evaluations
end

"""
    function gnewton(f::Function, df::Function, x::V, eps::Number; trace=false::Bool, itmax=10000::Int)

		Implements guarded Newton method for solving `f(x) == 0`,
### Arguments
-  `f(x)` is the function value at `x`,
-  `df(x)` is the Jacobian matrix of `f` at `x`,
-  `x` (input) is the initial guess,
-  `eps > 0` is the error tolerance: the method stops when `||f(x)|| < eps`.

### Optional named arguments
-  `xlist` is the list of the iterates (if `xlist != nothing`)
-  `trace=1` to trace the execuation of the algorithm
-  `itmax=m` to limit the algorithm to no more than `m` iterations (default is 10000)

  Returns solution, along with number of function and derivative (i.e., Jacobian) evaluations (`nfe`, `nde` resp.)

  Example: To solve ``x-\\cos(x) = 0`` with a tolerance of ``10^{-10}`` and starting at ``x=0``:
    ```x,nfe,nde = gnewton(x->x-cos(x),x->1+sin(x),0,1e-10)```
"""
function gnewton(f::Function, df::Function, x::V, eps::Number; 
        trace::Int=0, itmax::Int=10000, xlist::Union{Vector{V},Nothing} = nothing) where {V}

c1 = one(eltype(V))/10 # sufficient decrease parameter
fx = f(x)
nfe = 1
nde = 0
iter = 0
while norm(fx) >= eps && iter <= itmax
  dfx = df(x)
  nde = nde + 1
  d = - (dfx \ fx)
  alpha = one(eltype(V))
  fx_new = f(x+alpha*d)
  nfe = nfe + 1
  while norm(fx_new,2) >= (1-c1*alpha)*norm(fx,2) && alpha*norm(d) >= 1e-12*(1+norm(x))
      alpha = alpha / 2
      fx_new = f(x+alpha*d)
      nfe = nfe + 1
  end
  if alpha*norm(d) < 1e-12*(1+norm(x))
    return x, nfe, nde
  end
  if trace > 0; println("gnewton: alpha = ",alpha) end
  x = x + alpha*d
  fx = fx_new
  if trace > 1
    println("newton: iter = ", iter, ", x = ", x,", fx = ",fx,", dfx (old) = ",dfx)
  elseif trace > 0
    println("newton: iter = ", iter, ", ||f(x)|| = ",norm(fx))
  end
  if xlist != nothing; push!(xlist,x) end
  iter = iter + 1
end
return x,nfe,nde
end # function


""" 
    x,nfe = regfalsi(f::Function,a::T,b::T,eps::T;
        max_iter::Int=100, trace::Int=0, xlist::Union{Vector{T},Nothing}=nothing)

    Uses the Regula Falsi method on interval ``[a,b]`` (or ``[b,a]`` if ``a > b``) to
    compute approximate solution `x` to `f(x) == 0`.

### Arguments
- `f(x)` function 
- `a` 1st endpoint of interval
- `b` 2nd endpoint of interval
- `eps` target accuracy
### Optional named arguments
- `max_iter` maximum number of iterations (default is 100)
- `trace` trace execution of algorithm (default is 0 for no tracing)
- `xlist` a history of the iterates if `xlist != nothing`
"""
function regfalsi(f::Function,a::T,b::T,eps::T; 
    max_iter::Int=100, trace::Int=0, xlist::Union{Vector{T},Nothing}=nothing) where {T}
fa = f(a); fb = f(b);
nfe = 2

if sign(fa)*sign(fb) > 0
    println("regfalsi: f(a0) and f(b0) are of the same sign.  Stop!")
    return
end

c = a - fa*(b-a)/(fb-fa)
it_count = 0;
if trace > 0
    println(" iter\t  a\t  b\t  c\t  f(c)\t b-c")
end
while b-c > eps && it_count < max_iter
    it_count = it_count + 1;
    fc = f(c);
    nfe = nfe + 1
    if trace > 0
        println("regfalsi: ",it_count,"\t", a,"\t", b,"\t", c,"\t", fc,"\t\t", b-c)
    end
    if ( sign(fb)*sign(fc) <= 0 )
        a = c
        fa = fc
    else
        b = c
        fb = fc
    end
    c = a - fa*(b-a)/(fb-fa)
    if xlist != nothing
        push!(xlist,c) 
    end
end
return (c,nfe)
end # function


end # module 
