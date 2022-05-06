# code for Runge-Kutta solution of differential complementarity problems of index one.

"""
    module RK

Module for Runge-Kutta methods represented by a Butcher tableau in the `RKmethod` structure.
"""
module RK

# general functions
export RKmethod, rk, rk_explicit, rk_implicit # , rkfunc, drkfunc
export expEuler, impEuler, heun, trap, rk4, radauIIA2, radauIIA3, dirkA, dirkCR

# require("newton.jl")

# using Newton

"""
    RKmethod{T}

Runge-Kutta method represented by Butcher tableau.
`T` is the floating point type. 
The fields are:
- `A` the main coefficients ``a_{ij}`` of the tableau; an ``s\\times s`` matrix
- `b` the coefficients ``b_i`` for the final result of the method
- `c` the coefficients ``c_i`` for the time; note that ``c_i = \\sum_{j=1}^s a_{ij}``

The Runge-Kutta equations are
``v_i = f(t_n+c_i h, x_n+h\\sum_{j=1}^s a_{ij}v_j) `` and the output is
``x_{n+1} = x_n +h\\sum_{j=1}^s b_j v_j``.
"""
struct RKmethod{T}
  A :: Array{T,2}
  b :: Array{T,1}
  c :: Array{T,1}
  explicit :: Bool
end

# various RK methods
"""
    expEuler([T::Type])

Explicit Euler method Butcher tableau
"""
expEuler(T::Type) = RKmethod(zeros(T,1,1),[one(T)],[zero(T)], true)
expEuler()        = expEuler(Float64)
"""
    impEuler([T::Type])

Implicit Euler method Butcher tableau
"""
impEuler(T::Type) = RKmethod( ones(T,1,1),[one(T)],[one(T)], false)
impEuler()        = impEuler(Float64)
"""
    heun([T::Type])

Heun's method (explicit) Butcher tableau
"""
heun(T::Type)     = RKmethod([zero(T) zero(T); one(T) zero(T)],[one(T)/2;one(T)/2],[zero(T);one(T)], true)
heun()            = heun(Float64)
"""
    trapRK([T::Type])

Implicit trapezoidal method Butcher tableau
"""
trapRK(T::Type)   = RKmethod([zero(T) zero(T); one(T)/2;one(T)/2],[one(T)/2;one(T)/2],[zero(T);one(T)], false)
trapRK()          = trap(Float64)
"""
    rk4([T::Type])

Traditional 4th order Runge-Kutta method Butcher tableau
"""
rk4(T::Type)      = begin 
    z = zero(T); six = T(6); half = T(1)/T(2)
    RKmethod([z z z z; half z z z; z half z z; z z T(1) 0],
                           [T(1); T(2); T(2); T(1)]/six, [z;half;half;T(1)], true)
end
rk4()             = rk4(Float64)
"""
    radauIIA2([T::Type])

Two-stage Radau IIA method Butcher tableau (3rd order)
"""
radauIIA2(T::Type) = RKmethod([T(5)/T(12) T(-1)/T(12); T(3)/T(4) T(1)/T(4)],[T(3)/T(4);T(1)/T(4)],[T(1)/T(3),T(1)], false)
radauIIA2()        = radauIIA2(Float64)
"""
    radauIIA3([T::Type])

Three-stage Radau IIA method Butcher tableau
"""
radauIIA3(T::Type) = RKmethod([(88-7sqrt(T(6)))/360 (296-169sqrt(T(6)))/1800 (-2+3sqrt(T(6)))/225;
                      (296+169sqrt(T(6)))/1800 (88+7sqrt(T(6)))/360 (-2-3sqrt(T(6)))/225;
                      (16-sqrt(T(6)))/36 (16+sqrt(T(6)))/36 T(1)/9],
                     [(16-sqrt(T(6)))/36; (16+sqrt(T(6)))/36; T(1)/9],
                     [(4-sqrt(T(6)))/10; (4+sqrt(T(6)))/10; T(1)], false)
radauIIA3()        = radauIIA3(Float64)
# DIRK methods of Alexander (dirkA), and Crouzeix & Raviart (dirkCR)
"""
    dirkA([T::Type])

Alexander's DIRK method Butcher tableau
"""
dirkA(T::Type) = begin
    theta = atan(T(2)^(-T(3)/T(2)))/T(3)
    alpha = cos(theta)/sqrt(T(2)) + 1
    tau2 = (1+alpha)/2
    b1 = -(6alpha^2-16alpha+1)/4
    b2 =  (6alpha^2-20alpha+5)/4
    z = zero(T)
    RKmethod([alpha z z; tau2-alpha alpha z; b1 b2 alpha],
            [b1;b2;alpha],[alpha;tau2;1], false)
end
dirkA() = dirkA(Float64)
"""
    dirkCR([T::Type])

Crouzeix and Raviart's DIRK method Butcher tableau
"""
dirkCR(T::Type) = begin
    gamma = cos(T(pi)/18)/sqrt(T(3))+T(1)/T(2)
    delta = 1/(6*(2gamma-1)^2)
    RKmethod([gamma 0 0; 1/2-gamma gamma 0; 2gamma 1-4gamma gamma],
              [delta; 1-2delta; delta],[gamma; 1/2; 1-gamma], false)
end
dirkCR() = dirkCR(Float64)

"""
    isexplicit(rkm::RKmethod{T}) 

Returns true if rkm.A is strictly upper triangular
"""
isexplicit(rkm::RKmethod{T}) where {T} = 
    all(iszero(rkm.A[i,j]) for j = 1:size(rkm.A,1) for i = 1:j)
"""
    isdirk(rkm::RKmethod{T})

Returns true if rkm.A is upper triangular with all diagonal entries non-zero
"""
isdirk(rkm::RKmethod{T}) where {T} = 
    all(iszero(rkm.A[i,j]) for j = 1:size(rkm.A,1) for i = 1:j-1) &&
    ! any(iszero(rkm.A[i,i]) for i = 1:size(rkm.A,1))

"""
    rk(f::Function,rkm::RKmethod{T},t0::T,x0::V,h::T,n::Int;
        eps::T=T(10.0)^-6, df::Union{Function,Nothing}, solver!::Function=iter_solver!) where {T,V}
Solves ``dx/dt=f(t,x)`` with initial condition ``x(t_0)=x_0\\in V`` using Runge-Kutta method with tableau `rkm`.

Note that `V` is an approximate vector space over `T` and ``x(t)\\in V`` for all ``t``.

### Arguments
- `x0` is the initial value
- `f(t,x)` is the function defining the differential equation
- `h` is the step size for the method
- `n` is the number of steps to take of the method
### Optional named arguments
- `eps`

Returns the pair `(xlist,tlist)`

"""
function rk(f::Function,rkm::RKmethod{T},t0::T,x0::V,h::T,n::Int;
        eps::T=T(10.0)^-6, df::Union{Function,Nothing}, solver!::Function=iter_solver!) where {T,V}
    

end

"""
    function iter_solver!(f::Function,df::Function,x::V,vlist::Vector{V},A::Array{T,2},eps::T;
    maxit::Int=1000) where {T,V}

This is a simple (fixed-point) iterative solver for the Runge-Kutta equations: 
for ``i=1,\\ldots,s``:
``v_i\\gets f(t_i+c_ih,x+\\sum_{j<i}a_{ij}v_j)``
This is, in effect, a block nonlinear Gauss-Seidel method. It repeats until the maximum change in ``v_i`` is no more than `eps` or `maxit` iterations are reached.
In the case of explicit methods, one iteration is sufficient to obtain the solution.
"""
function iter_solver!(f::Function,df::Function,x::V,vlist::Vector{V},A::Array{T,2},eps::T;
    maxit::Int=1000) where {T,V}
    
    firstloop = true; maxerr = zero(T)
    A = rkm.A; c = rkm.c
    while firstloop || maxerr > eps
        maxerr = zero(T)
        for i = 1:size(A,1)
            xin = x
            for j = 1:size(A,2)
                xin += (h*A[i,j])*vlist[j]
            end
            temp = f(t+c[i]*h,xin)
            errval = norm(temp-vlist[i])
            maxerr = max(maxerr,erval)
            vlist[i] = temp
        end
        firstloop = false
    end
    vlist
end

"""
    rk_explicit(f::Function,rkm::RKmethod{T},t0::T,x0::V,h::T,n::Int) where {T,V}

Solves ``dx/dt=f(t,x)`` with initial condition ``x(t_0)=x_0\\in V`` using an explicit Runge-Kutta method with tableau `rkm`.

Note that `V` is an approximate vector space over `T` and ``x(t)\\in V`` for all ``t``.

### Arguments
- `x0` is the initial value
- `f(t,x)` is the function defining the differential equation
- `h` is the step size for the method
- `n` is the number of steps to take of the method

Returns the pair `(xlist,tlist)`
"""
function rk_explicit(f::Function,rkm::RKmethod{T},t0::T,
                        x0::V,h::T,n::Int) where {T,V}
  # Assume that the method is explicit and rk.A is strictly lower triangular
  # println("RK method: ",rkm,", t0 = ",t0,", x0 = ",x0,", h = ",h,", n = ",n)
  x = x0
  t = t0
  xlist = V[]
  tlist = T[]
  push!(xlist,x)
  push!(tlist,t)
  v  = Vector{V}(undef,size(rkm.A,1))
  fv = Vector{V}(undef,size(rkm.A,1))
  s = size(rkm.A,2) # number of stages
  for step = 1:n
    for i = 1:s
      if i > 1
        fv[i-1] = f(t+rkm.c[i-1]*h,v[i-1])
      end
      v[i] = x
      for j = 1:i-1
        v[i] += h*rkm.A[i,j]*fv[j]
      end
    end
    fv[s] = f(t+rkm.c[s]*h,v[s])
    new_x = x
    for j = 1:s
      new_x += h*rkm.b[j]*fv[j]
    end
    # update xlist, tlist, x, t, etc.
    x = new_x
    t = t + h
    push!(xlist,x)
    push!(tlist,tlist[end] + h)
  end
  return xlist,tlist
end

# Runge-Kutta method for solving dx/dt = f(t,x), x(t0) = x0
# takes n steps of length h, solving the RK equations to an accuracy of eps or better
# returns (tlist,xlist)

function rk_implicit(f::Function,df::Function,rkm::RKmethod{T},t0::T,
                        x0::V,h::T,n::Int,eps::T; trace=false::Bool) where {T,V}
  # Assume that the method is explicit and rk.A is strictly lower triangular
  # println("RK method: ",rkm,", t0 = ",t0,", x0 = ",x0,", h = ",h,", n = ",n)
  x = x0
  t = t0
  xlist = Array(T,length(x),n+1)
  tlist = Array(T,n+1)
  xlist[:,1] = x
  tlist[1] = t
  v  = Array(T,length(x),size(rkm.A,1))
  F(vvec) = rkfunc(f,rkm,t,x,h,vvec)
  dF(vvec) = drkfunc(f,df,rkm,t,x,h,vvec)
  s = size(rkm.A,2) # number of stages
  for step = 1:n
    vvec,nFe,ndFe = newton(vvec->rkfunc(f,rkm,t,x,h,vvec),
                vvec->drkfunc(f,df,rkm,t,x,h,vvec),v[:],eps)
    v = reshape(vvec,length(x0),s)
    new_x = x
    for j = 1:s
      new_x += h*rkm.b[j]*f(t+rkm.c[j]*h,v[:,j])
    end
    # update xlist, tlist, x, t, etc.
    x = new_x
    t = t + h
    xlist[:,step+1] = x
    tlist[step+1] = tlist[step] + h
  end
  return tlist,xlist
end

"""
    rkfunc(f::Function,rkm::RKmethod{T},t::T,x::V,h::T,vvec::Array{V,1}) where {T,V}
Function giving Runge-Kutta equations to solve (i.e., set output to zero)
"""
function rkfunc(f::Function,rkm::RKmethod{T},t::T,
            x::V,h::T,vvec::Array{V,1}) where {T,V}
  # v = reshape(vvec,length(x),div(length(vvec),length(x)))
  out = Array{V}(undef,length(vvec))
  for i = 1:size(rkm.A,1)
    out[i] = v[i] - x
    for j = 1:size(rkm.A,2)
      out[i] -= h*rkm.A[i,j]*f(t+h*rkm.c[j],v[j])
            
    end
  end
  return out
end

# Function returning Jacobian of R-K equations to solve
function drkfunc(f::Function,df::Function,rkm::RKmethod{T},t::T,
            x::Array{T,1},h::T,vvec::Array{T,1}) where {T}
  v = reshape(vvec,length(x),div(length(vvec),length(x)))
  Jout = zeros(T,length(vvec),length(vvec))
  s = size(rkm.A,1)
  for i = 1:s
    iblk = [(i-1)*length(x)+1:i*length(x)]
    Jout[iblk,iblk] = eye(length(x))
    for j = 1:s
      jblk = [(j-1)*length(x)+1:j*length(x)]
      Jout[iblk,jblk] -= h*rkm.A[i,j]*df(t+h*rkm.c[j],v[:,j])
    end
  end
  return Jout
end
end # module RK
