"""
	module CSpline
	Cubic spline interpolation module
"""
module CSpline

export CSpl, ev, evd, evdd, cspline_nat, cspline_clamped, cspline_nak
using LinearAlgebra

"""
	struct CSpl{T}

	Fields:
	x::Array{T,1} # x[i] is i'th breakpoint
	M::Array{T,1} # 2nd derivative values at breakpoints
	C::Array{T,1}
	D::Array{T,1}

On [x[i],x[i+1]], the spline function is
	s(t) = M[i]*(x[i+1]-t)^3/(6h[i])+M[i+1]*(t-x[i])^3/(6h[i])
			+ C[i]*(x[i+1]-t) + D[i]*(t-x[i])
where h[i] = x[i+1]-x[i].
"""
struct CSpl{T,V}
	x::Array{T,1} # x[i] is i'th breakpoint
	M::Array{V,1} # 2nd derivative values at breakpoints
	C::Array{V,1}
	D::Array{V,1}
end # type

"""
	ev(s,t) evaluates spline s at t: s(t)
"""
function ev(s::CSpl{T,V}, t::T)::V where {T,V}
	n = length(s.x)
	i = searchsortedlast(s.x,t)
	if i < 1
		i = 1
	end
	if i >= n
		i = n-1
	end
	z1 = s.x[i+1] - t; z2 = t - s.x[i]; h = s.x[i+1]-s.x[i]
	y = (s.M[i]*z1^3 + s.M[i+1]*z2^3)/(6h) + s.C[i]*z1 + s.D[i]*z2
end # function

"""
	evd(s,t) evaluates spline derivative s' at t: s'(t)
"""
function evd(s::CSpl{T,V}, t::T)::V where {T,V}
	n = length(s.x)
	i = searchsortedlast(s.x,t)
	if i < 1
		i = 1
	end
	if i >= n
		i = n-1
	end
	z1 = s.x[i+1] - t; z2 = t - s.x[i]; h = s.x[i+1]-s.x[i]
	y = (-s.M[i]*z1^2 + s.M[i+1]*z2^2)/(2h) - s.C[i] + s.D[i]
end # function

"""
	evdd(s,t) evaluates spline 2nd derivative s'' at t: s''(t)
"""
function evdd(s::CSpl{T,V}, t::T)::V where {T,V}
	n = length(s.x)
	i = searchsortedlast(s.x,t)
	if i < 1
		i = 1
	end
	if i >= n
		i = n-1
	end
	z1 = s.x[i+1] - t; z2 = t - s.x[i]; h = s.x[i+1]-s.x[i]
	y = (s.M[i]*z1 + s.M[i+1]*z2)/h
end # function

# Array versions
ev(s::CSpl{T,V}, ta::Array{T,N}) where {T,V,N} = (t->ev(s,t)).(ta)
evd(s::CSpl{T,V}, ta::Array{T,N}) where {T,V,N} = (t->evd(s,t)).(ta)
evdd(s::CSpl{T,V}, ta::Array{T,N}) where {T,V,N} = (t->evdd(s,t)).(ta)

include("tridiag.jl")

import .TDG: Tridiag, lufact!, lusolve!

function order_xy!(x::Array{T,1}, y::Array{V,1}) where {T,V}
	# put the x-array into order, and permute y consistently
	p = sortperm(x)
	x = x[p]
	y = y[p]
end

function base_system(x::Array{T,1},y::Array{V,1}) where {T,V}
	# creates base tridiagonal linear system
	if length(x) != length(y)
		raise(DimensionMismatchError())
	end
	n = length(x)
	a1 = (x[3:n]-x[1:n-2])/3
	b1 = (x[2:n]-x[1:n-1])/6
	c1 = (x[2:n]-x[1:n-1])/6
	td = Tridiag([0;a1;0],b1,c1)
	rhs = [zero(V);(y[3:n]-y[2:n-1])./(x[3:n]-x[2:n-1]) - (y[2:n-1]-y[1:n-2])./(x[2:n-1]-x[1:n-2]);zero(V)]
	return td,rhs
end # function

function compute_CD(x::Array{T,1},y::Array{V,1},M::Array{V,1}) where {T,V}
	# computes C, D arrays once we have computed M from x & y
	C = zeros(V,length(x))
	D = zeros(V,length(x))
	for i = 1:length(x)-1
		h = x[i+1] - x[i]
		C[i] = (y[i  ]/h) - M[i  ]*h/6
		D[i] = (y[i+1]/h) - M[i+1]*h/6
	end
	return C, D
end # function

"""
	cs = cspline_nat(x,y)
Computes natural cubic spline s where s(x[i]) = y[i], i = 1, ..., length(x),
and s''(x[1]) == s''(x[length(x)]) == 0.
"""
function cspline_nat(x::Array{T,1},y::Array{V,1}) where {T,V}
	if length(x) == 0
		throw(ArgumentError("x[] is an empty array"))
	end
	if length(x) != length(y)
		throw(DimensionMismatch("length(x) != length(y)"))
	end
	x = copy(x); y = copy(y)
	order_xy!(x,y)
	td, M = base_system(x,y)
	t = x[1] # type of t is T
	td.a[1]   = one(t)
	td.a[end] = one(t)
	td.b[end] = zero(t)
	td.c[1]   = zero(t)
	M[1]   = zero(t)
	M[end] = zero(t)
	lufact!(td)
	lusolve!(td,M)
	C, D = compute_CD(x,y,M)
	return CSpl(x,M,C,D)
end # function

"""
	cs = cspline_clamped(x,y,y1p,y2p)
Computes clamped spline s where s(x[i]) = y[i], i = 1, ..., length(x),
and s'(x[1]) == y1p, and s'(x[length(x)]) == y2p.
"""
function cspline_clamped(x::Array{T,1},y::Array{V,1},y1p::V,y2p::V) where {T,V}
	if length(x) == 0
		throw(ArgumentError("x[] is an empty array"))
	end
	if length(x) != length(y)
		throw(DimensionMismatch("length(x) != length(y)"))
	end
	x = copy(x); y = copy(y)
	order_xy!(x,y)
	td, M = base_system(x,y)
	n = length(x)
	h1 = x[2]-x[1]
	td.a[1] = h1/3
	td.c[1] = h1/6
	M[1] = (y[2]-y[1])/h1 - y1p
	h2 = x[n]-x[n-1]
	td.a[n] = h2/3
	td.b[n-1] = h2/6
	M[n] = y2p - (y[n]-y[n-1])/h2
	lufact!(td)
	lusolve!(td,M)
	C, D = compute_CD(x,y,M)
	return CSpl(x,M,C,D)
end # function

"""
	cs = cspline_nak(x,y)
Computes "not-a-knot" spline s where s(x[i]) = y[i], i = 1, ..., length(x),
and s''' continuous at x[2] and x[length(x)-1].
"""
function cspline_nak(x::Array{T,1},y::Array{V,1}) where {T,V}
	if length(x) == 0
		throw(ArgumentError("x[] is an empty array"))
	end
	if length(x) != length(y)
		throw(DimensionMismatch("length(x) != length(y)"))
	end
	x = copy(x); y = copy(y)
	order_xy!(x,y)
	td, M = base_system(x,y)
	n = length(x)
	# We need to use the Sherman-Morrisson-Woodbury formula
	# as we need to solve a rank-2 modification of a tridiagonal system.
	h = x[2:n] - x[1:n-1]
	td.a[1] = h[2]
	td.c[1] = -(h[1]+h[2])
	td.a[n] = h[n-2]
	td.b[n-1] = -(h[n-1]+h[n-2])
	# Matrix is td + [h[1]*e_1, h[n-2]*e_n]*[e_3, e_{n-2}]' = td + U*V'
	# Inverse of this matrix is
	# inv(td+U*V') = inv(td) - inv(td)*U*inv(eye+V'*inv(td)*U)*V'*inv(td)
	t = x[1] # used for type information later
	lufact!(td) # factorize td
	U = zeros(T,n,2); W = zeros(T,n,2)
	U[1,1] = h[1]
	U[n,2] = h[n-1]
	W[3,1] = one(T)
	W[n-2,2] = one(T)
	lusolve!(td,U) # now U is solution of tridiag system
	lusolve!(td,M) # now M is solution of tridiag system
	M = M - U*((W'*U+I)\(W'*M))
	C, D = compute_CD(x,y,M)
	return CSpl(x,M,C,D)
end

end # module
