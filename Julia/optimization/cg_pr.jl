"""
 [x,fval,gval,nfe,nge,xs] = cg_pr(func,dfunc,x,tol,c1,c2,it_limit;trace=0)
 Polak-Ribiere conjugate gradient method (without restart)
 Uses Wolfe line-search conditions (parameters c1 and c2)

 func is the name of the function to minimize: func(x)
 dfunc is the name of the gradient function: dfunc(x)
 x is the current position
 tol is the stopping tolerance: ||gradient|| < tol
 c1 and c2 are the Wolfe line-search parameters: 0 < c1 < c2 < 1
 it_limit is the iteration limit
 if trace is non-zero, then print out information about process

 Returns
	x	-- final "minimizer"
	fval	-- func(x), function value at x
	gval	-- gradient value at x
	nfe	-- # function eval'ns
	nge	-- # gradient eval'ns
	xs	-- list of iterates (only used if trace ~= 0)
"""
function cg_pr(func::Function,dfunc::Function,x::V,tol::T,c1::T,c2::T;
        it_limit::Integer=10000,trace::Integer=0) where {T,V}

xs = V[];
n = length(x);
alpha0 = one(T)

fval = func(x);
g    = dfunc(x)
nfe = 1;
nge = 1;
p = -g;
g_sqr = g'*g;
for k = 0:it_limit
  if  trace > 0
    # xs = [xs;[x]];
    push!(xs,x)
  end
  if trace > 1
    println("cg_pr: cos(theta) = ", -dot(p,g)/(norm(p)*norm(g)))
    println("cg_pr: ||p||/||g|| = ", norm(p)/norm(g))
  end
  println("Starting wolve_search2")

  alpha,fval,new_g,nfe_ws,nge_ws =
        wolfe_search2(func,dfunc,x,p,alpha0,c1,c2,fval,g;trace=max(trace-1,0));
  nfe = nfe + nfe_ws;
  nge = nge + nge_ws;
  if trace > 0
    println("cg_pr: iter # ",k,", alpha = ",alpha,", function value = ", fval);
  end
  x = x + alpha*p;
  new_g_sqr = new_g'*new_g;
  # Polak Ribiere modification
  beta = new_g'*(new_g-g)/g_sqr;
  g_sqr = new_g_sqr;
  g = new_g;
  if norm(g,2) < tol
    return x,fval,g,nfe,nge,xs
  end
  # Restart code (if desired)
  # if k # n == n-1
  #   beta = 0;
  # end
  p = -g + beta*p;
end
println("Have used ",it_limit," iterations; returning without finding solution");
return x,fval,g,nfe,nge,xs
end # function
