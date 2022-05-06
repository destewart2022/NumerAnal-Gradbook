"""
    wolfe_search2(f,df,x,p,alpha0,c1,c2,f0,df0,trace) --> alpha,fval,dfval,nfe,nge
Performs search for point satisfying the two strong Wolfe conditions:

f(x+alpha p) <= f(x) + c1.alpha p'.grad f(0)		(WC1)
|p'.grad f(x+alpha p)| <= c2|p'.grad f(x)|	    	(WC2b)

alpha0 is the initial value of alpha in the search.  If alpha0
	satisfies the Wolfe conditions, then the search stops there.
f0 is f(x) at the initial x.
df0 is grad f(x) at the initial x.
f is the objective function.
df is the gradient function.
If trace is non-zero then print out information about the process

The algorithm follows that of Wright & Nocedal "Numerical Optimization", pp. 60-62, section 3.4

Returns
	alpha	-- step length parameter
	fval	-- function value at end
	dfval	-- gradient value at end
	nfe     -- # function evaluations
	nge	    -- # gradient evaluations

Note: V is an inner product space over R. Then f:V->R and df:V->V.
"""
function wolfe_search2(f::Function,df::Function,x::V,p::V,alpha0::R,
    c1::R,c2::R,f0::R,df0::V;trace::Integer=0) where {R,V}

if trace > 0
    println("Entered wolfe_search2")
end
#
# Check parameters
#
alpha = 0;
fval = f0;
dfval = df0;
nfe = 0;
nge = 0;

if ! ( zero(R) < c1 < c2 < one(R) )
  println("wolfe_search2: Error: Need 0 < c1 < c2 < 1");
  if trace > 0
    println("wolfe_search2: c1 = $c1, c2 = $c2");
  end
  return alpha,fval,dfval,nfe,nge
end

# Initialization
slope0 = dot(p,df0);
if slope0 >= 0
  println("wolfe_search2: Error: Need a descent direction");
  if trace > 0
    println("wolfe_search2: p'.grad f(x) = ",slope0);
  end
  return alpha,fval,dfval,nfe,nge
end

#
# Bracketing phase
#
if trace > 0
  println("wolfe_search2: Bracketing phase: alpha = ", alpha0);
end
alpha = alpha0;
old_alpha = 0;
old_fval  = f0;
old_dfval = df0;
old_slope = slope0;
if trace > 0
  println("wolfe_search2: f(x) = $f0, p'.grad f(x) = ",slope0);
end

# Main loop
firsttrip = true;
while true ### forever do...
  xplus = x+alpha*p;
  fval = f(xplus)
  dfval = df(xplus)
  nfe = nfe + 1;
  nge = nge + 1;
  if trace > 0
    println("wolfe_search2: alpha = $alpha, f(x+alpha*p) = ", fval);
  end
  if ( fval > f0 + c1*alpha*slope0 ) || ( ( ! firsttrip ) &&
	( fval >= old_fval ) )
    if trace > 0
      println("wolfe_search2: (WC1) failed or f increased");
    end
    break;
  end
  if trace > 0
    println("wolfe_search2: (WC1) holds & f decreased");
  end
  slope = dot(p,dfval);
  if trace > 0
    println("wolfe_search2: p'.grad f(x+alpha*p) = ", slope);
  end
  if ( abs(slope) <= c2*abs(slope0) )
    if trace > 0
      println("wolfe_search2: (WC2) holds");
    end
    return alpha,fval,dfval,nfe,nge
  end
  if ( slope >= 0 )
    if trace > 0
      println("wolfe_search2: f'(alpha) >= 0");
    end
    break;
  end

  # Update variables -- note no upper limit on alpha
  temp = alpha;
  # alpha = alpha + old_alpha;
  alpha = 2*alpha;
  old_alpha = temp;
  old_fval = fval;
  old_slope = slope;
end

if ( trace > 0 )
  println("wolfe_search2: Entering ''zoom'' phase.")
end

#
# "Zoom" phase
#
alpha_lo = old_alpha;
alpha_hi =     alpha;
f_lo     = old_fval;
f_hi     =     fval;
df_lo    = old_dfval;
slope_lo = old_slope;
#df_hi   =     dfval;
#slope_hi=     slope;

if trace > 0
  println("wolfe_search2: zoom phase: alpha_lo = $alpha_lo, alpha_hi = %$alpha_hi");
end

iter_cnt = 0;
while abs(alpha_hi-alpha_lo) > 1e-15*alpha0

  # form quadratic interpolant of function values on alpha_lo, alpha_hi
  # and the derivative at alpha_lo...
  # and find min of interpolant within interval [alpha_lo,alpha_hi]
  a = f_lo;
  b = slope_lo;
  dalpha = alpha_hi-alpha_lo;
  c = (f_hi - f_lo - dalpha*slope_lo)/dalpha^2;
  # iter_cnt
  if ( ( c <= 0 ) | ( mod(iter_cnt,3) == 2 ) )
    # Use bisection
    alpha = alpha_lo + 0.5*dalpha;
    if trace > 0
      println("wolfe_search2: using bisection: c=", c);
    end
  else
    # Use min of quadratic
    alpha = alpha_lo - 0.5*b/c;
    if trace > 0
      println("wolfe_search2: using quadratic: a=$a, b=$b, c=$c");
    end
  end
  if trace > 0
    println("wolfe_search2: alpha = ", alpha);
  end

  # main part of loop
  xplus = x + alpha*p;
  fval = f(xplus);
  dfval = df(xplus);
  nfe = nfe + 1;
  nge = nge + 1;
  if trace > 0
    println("wolfe_search2: f(x+alpha*p) = ", fval);
  end
  if ( ( fval > f0 + c1*alpha*slope0 ) | ( fval >= f_lo ) )
    if trace > 0
      println("wolfe_search2: zoom: (WC1) fails or f increased");
      println("wolfe_search2: zoom: update alpha_hi");
    end
    alpha_hi = alpha;
    f_hi   = fval;
    # df_hi = feval(df,xplus);
  else
    if trace > 0
      println("wolfe_search2: zoom: (WC1) holds and f decreased");
    end
    slope = dot(p,dfval);
    if trace > 0
      println("wolfe_search2: p'.grad f(x+alpha*p) = ", slope);
    end
    if ( abs(slope) <= c2*abs(slope0) )
      if trace > 0
	    println("wolfe_search2: zoom: (WC2b) holds");
      end
      return alpha,fval,dfval,nfe,nge
    end
    if ( slope*dalpha >= 0 )
      if trace > 0
	    println("wolfe_search2: zoom: alpha_hi <- alpha_lo & update alpha_lo");
      end
      alpha_hi = alpha_lo;
      f_hi     = f_lo;
      # df_hi  = df_lo;
      # slope_hi= slope_lo;
      alpha_lo = alpha;
      f_lo     = fval;
      df_lo    = dfval;
      slope_lo = slope;
    else
      if trace > 0
	    println("wolfe_search2: zoom: update alpha_lo");
      end
      alpha_lo = alpha;
      f_lo     = fval;
      df_lo    = dfval;
      slope_lo = slope;
    end
  end
  # Update iteration count
  iter_cnt = iter_cnt + 1;

end # of forever do...

if trace > 0
    println("Leaving wolfe_search2")
end

return alpha,fval,dfval,nfe,nge
end # function
