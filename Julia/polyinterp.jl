module PolyInterp

export divdif, dd_interp, dd_hermite, dd_hinterp

"""
	dd = divdif(xs::Array{T,1},ys::Array{V,1})
Construct divided difference table

xs = array of interpolation points;
ys = array of (function) values to interpolate

Returns dd, the divided difference table
See also dd_interp().
Note that V must be a (approximate) vector space over an
approximate field T.
"""
function divdif(xs::Array{T,1},ys::Array{V,1}) where {T,V}
  n = length(xs) - 1;
  dd = copy(ys)

  for k = 1:n
    # compute k'th order divided differences
    for j = n:-1:k
      dd[j+1] = (dd[j+1] - dd[j])/(xs[j+1] - xs[j-k+1]);
    end
  end
  dd
end # function

"""
	dd = dd_hermite(xs::Array{T,1},fs::Array{V,1},dfs::Array{V,1})
Computes divided difference table for Hermite interpolation:
interpolating polynomial satisfies

p(xs[i]) = fs[i],
p'(xs[i]) = dfs[i]

See also dd_hinterp.m.
Note that V must be a (approximate) vector space over an
approximate field T.
"""
function dd_hermite(xs::Array{T,1},fs::Array{V,1},dfs::Array{V,1}) where {T,V}
  n = length(xs) - 1;
  dd = zeros(2*n+2);
  dd[1] = fs[1];
  dd[2] = dfs[1];
  for i = 1:n
    dd[2*i+1] = (fs[i+1]-fs[i])/(xs[i+1]-xs[i]);
    dd[2*i+2] = dfs[i+1];
  end

  new_xs = zeros(1,2*length(xs));
  new_xs[2*(1:length(xs))]   = xs;
  new_xs[(2*(1:length(xs))).-1] = xs;
  for i = 2:(2*n+1)
    for j = (2*n+1):-1:i
      dd[j+1] = (dd[j+1] - dd[j]) / (new_xs[j+1] - new_xs[j-i+1]);
    end
  end
  dd
end # function

"""
	pval  = dd_interp(xs::Array{T,1},dd::Array{V,1},t::U)
	pvals = dd_interp(xs::Array{T,1},dd::Array{V,1},t::Array{U,N})

Evaluates interpolations polynomial defined
by interpolation points xs and the divided difference table dd.

Evaluate the interpolating polynomial at t

See also divdif().
Note that V must be a (approximate) vector space over an
approximate field T, while U is T or a ring containing T.
"""
function dd_interp(xs::Array{T,1},dd::Array{V,1},t::U) where {T,U,V}
  n = length(xs) - 1;
  pval = dd[n+1];
  for i = n-1:-1:0
    pval = dd[i+1] + (t-xs[i+1])*pval;
  end
  pval
end # function

function dd_interp(xs::Array{T,1},dd::Array{V,1},t::Array{U,N}) where {T,U,N,V}
  # n = length(xs) - 1;
  # pvals = ones(T,size(t))*dd[n+1];
  # for i = n-1:-1:0
  #   pvals = dd[i+1] .+ (t.-xs[i+1]).*pvals;
  # end
  # pvals
  [dd_interp(xs,dd,tv) for tv in t]
end # function

"""
	pval  = dd_hinterp(xs::Array{T,1},dd::Array{V,1},t::U)
	pvals = dd_hinterp(xs::Array{T,1},dd::Array{V,1},t::Array{U,N})
Evaluates interpolations polynomial defined
by interpolation points xs and
the divided difference table dd
Evaluate the interpolating polynomial at t
If t is a vector or matrix, then the interpolating
polynomial is applied componentwise.
See also dd_interp.m and dd_hermite.m.
"""
function dd_hinterp(xs::Array{T,1},dd::Array{V,1},t::U) where {T,U,V}
  new_xs = zeros(T,2*length(xs));
  new_xs[2*(1:length(xs))]   = xs;
  new_xs[2*(1:length(xs)).-1] = xs;
  pval = dd_interp(new_xs,dd,t);
  pval
end # function

function dd_hinterp(xs::Array{T,1},dd::Array{V,1},t::Array{U,N}) where {T,U,N,V}
  new_xs = zeros(2*length(xs));
  new_xs[2*(1:length(xs))]    = xs;
  new_xs[2*(1:length(xs)).-1] = xs;
  dd_interp(new_xs,dd,t);
  # [dd_hinterp(xs,dd,tv) for tv in t]
end # function

end # module PolyInterp
