"""
    inv_iter(A,x0,mu0,num_iter) --> (lambdalist,xlist)

Performs inverse iteration to find an eigenvalue of A
* `x0`	-- starting eigenvector estimate
* `mu0`	-- starting eigenvalue estimate
* `num_iter` -- number of iterations
Output:
* `lambdalist`	-- list of eigenvalue estimates
* `xlist`		-- list of eigenvector estimates
"""
function inv_iter(A,x0,mu0,num_iter::Int)

    n = length(x0);

    x = x0;
    mu = mu0;
    lambdalist = zeros(typeof(mu0),num_iter+1);
    xlist      = Array{typeof(x0)}(undef, num_iter+1);
    lambdalist[1] = mu;
    xlist[1]   = x;

    for k = 1:num_iter
      y = (A-mu*I) \ x;
      yval,idx = findmax(abs.(x));
      # Estimate eigenvalue of (A-mu*I)^{-1}
      sigma = y[idx] / x[idx]
      mu = mu + 1/sigma;
      x  = y / norm(y);
      lambdalist[k+1] = mu;
      xlist[k+1]      = x;
    end
    return lambdalist,xlist
end # function

function zero(x::Array{T,N}) where {T,N}
    return fill(zero(T),size(x))
end