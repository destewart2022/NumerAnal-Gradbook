"""
	module TDG
	Includes Tridiag type
"""
module TDG

using LinearAlgebra
import Base.*
# import LinearAlgebra.lu, LinearAlgebra.lufact!

export Tridiag, lufact!, lusolve!

"""
		struct Tridiag{T}
		tridiagonal matrix with entries of type T
"""
struct Tridiag{T}
	# (length of b) = (length of c) = (length of a) - 1
	a::Array{T,1} # main diagonal
	b::Array{T,1} # sub-diagonal
  	c::Array{T,1} # super-diagonal
end # struct

Tridiag(n::Integer,T) = Tridiag(zeros(T,n),zeros(T,n-1),zeros(T,n-1))
Tridiag(n::Integer) = Tridiag(n,Float64)

"""
	lufact!(td::Tridiag)
Returns LU factorization in compact form; overwrites td.
No pivoting
"""
function lufact!(td::Tridiag{T}) where {T}
	n = length(td.a)
	if td.a[1] == 0
		throw(OverflowError())
	end
	td.b[1] = td.b[1] / td.a[1];
	for i = 2:n-1
      td.a[i] = td.a[i] - td.b[i-1]*td.c[i-1];
      if td.a[i] == 0
      	throw(OverflowError())
      end
      td.b[i] = td.b[i] / td.a[i];
    end
    td.a[n] = td.a[n] - td.b[n-1]*td.c[n-1];
    td
end # lu

function *(td::Tridiag{T},x::Array{T,1}) where {T}
	n = length(x)
	# println("n = ", n)
	if length(td.a) != length(x)
		throw(DimensionMissmatch())
	end
	y = td.a .* x + [zero(T);td.b .* x[1:n-1]] + [td.c .* x[2:n];zero(T)]
end # function

function *(td::Tridiag{T},x::Array{T,2}) where {T}
	n = size(x,1)
	# println("n = ", n)
	if length(td.a) != n
		throw(DimensionMissmatch())
	end
	y = zeros(eltype(x),size(x))
	for j = 1:size(x,2)
		y[:,j] = td.a .* x[:,j] + [zero(T);td.b .* x[1:n-1,j]] +
			 [td.c .* x[2:n,j];zero(T)]
	end
	y
end # function


"""
	lusolve!{T,V}(td::Tridiag{T}, y::Array{V,1})
	lusolve!{T}(td::Tridiag{T}, y::Array{T,2})

Solve tridiagonal system td*x = y (returning x) where td has been factored by lu(::Tridiag).
"""
function lusolve!(td::Tridiag{T}, y::Array{V,1}) where {T,V}
	n = length(td.a)
	if length(y) != n
		throw(DimensonMismatch())
	end
	for i = 2:n
		y[i] = y[i] - td.b[i-1]*y[i-1]
	end
	y[n] = y[n] / td.a[n]
	for i = n-1:-1:1
		if td.a[i] == 0
			throw(OverflowError())
		end
		y[i] = (y[i] - td.c[i]*y[i+1]) / td.a[i]
	end
	y
end # function

function lusolve!(td::Tridiag{T}, y::Array{T,2}) where {T}
	n = length(td.a)
	if size(y,1) != n
		throw(DimensonMismatch())
	end
	for i = 2:n
		y[i,:] = y[i,:] - td.b[i-1]*y[i-1,:]
	end
	y[n,:] = y[n,:] / td.a[n]
	for i = n-1:-1:1
		if td.a[i] == 0
			throw(OverflowError())
		end
		y[i,:] = (y[i,:] - td.c[i]*y[i+1,:]) / td.a[i]
	end
	y
end # function

end # module
