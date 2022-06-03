"""
    module Integration

A collection of integration methods.
"""
module Integration

export rectangle, trapezoidal, simpson, romberg

"""
    rectangle(f::Function,a::T,b::T,n::Int) where {T}

Uses rectangle rule with ``n`` function evaluations to estimate ``\\int_a^b f(x)\\,dx``
"""
function rectangle(f::Function,a::T,b::T,n::Int) where {T}
    h = (b-a)/n
    h*sum(f(a+i*h) for i=0:n-1)
end

"""
    trapezoidal(f::Function,a::T,b::T,n::Int) where {T}

Uses trapezoidal rule with ``n+1`` function evaluations to estimate ``\\int_a^b f(x)\\,dx``
"""
function trapezoidal(f::Function,a::T,b::T,n::Int) where {T}
    h = (b-a)/n
    h*((f(a)+f(b))/2 + sum(f(a+i*h) for i = 1:n-1))
end

"""
    simpson(f::Function,a::T,b::T,n::Int) where {T}

Uses Simpson's rule with ``2n+1`` function evaluations to estimate ``\\int_a^b f(x)\\,dx``
"""
function simpson(f::Function,a::T,b::T,n::Int) where {T}
    h = (b-a)/n
    (h/6)*(f(a)+f(b)+4*sum(f(a+(2i+1)*(h/2)) for i = 0:n-1) +
        2*sum(f(a+i*h) for i = 1:n-1))
end

"""
    gen_int_rule(f::Function,a::T,b::T,wts,nodes[,n::Int]) where {T}
    gen_int_rule(f::Function,a::T,b::T,wtnodes::Tuple{Any,Any}[,n::Int]) where {T}
    
Returns ``(b-a)\\sum_{i=1}^n w_i\\,f(a+(b-a)x_i)`` to approximate ``\\int_a^b f(t)\\,dt``; assumes that the weights ``w_i=`` `wts[i]` and nodes ``x_i=`` `nodes[i]` are for the interval ``[0,1]``.
    
    If `n` is given, then the composite version of the rule is used with `n` pieces, each of length ``(b-a)/n``.

    If the second form is used, then `wtnodes` is the pair `(wts,nodes)`.
"""
function gen_int_rule(f::Function,a::T,b::T,wts,nodes) where {T}
    sum(w*f(a+x*(b-a)) for (w,x) in zip(wts,nodes))*(b-a)
end
    
function gen_int_rule(f::Function,a::T,b::T,wts,nodes,n::Int) where {T}
    s = zero(f(a))
    for i =0:n-1
        s += sum(w*f(a+(i+x)*(b-a)/n) for (w,x) in zip(wts,nodes))*(b-a)/n
    end
    return s
end

function gen_int_rule(f::Function,a::T,b::T,wtnodes::Tuple{Any,Any}) where {T}
    gen_int_rule(f,a,b,wtnodes[1],wtnodes[2])
end
    
function gen_int_rule(f::Function,a::T,b::T,wtnodes::Tuple{Any,Any},n::Int) where {T}
    gen_int_rule(f,a,b,wtnodes[1],wtnodes[2],n)
end

using LinearAlgebra

"""
    romberg(f::Function,a::T,b::T,nlevels::Int) where {T}

Uses Romberg's method of estimating the integral of `f` over `[a,b]`.
The number of levels used is `nlevels` and the number of function evaluations is `2^{nlevals+1}-2`
"""
function romberg(f::Function,a::T,b::T,nlevels::Int) where {T}
    n = nlevels
    rtable = zeros(n,n);

    # Evaluate trapezoidal integrals

    # Evaluate simple trapezoidal rule
    rtable[1,1] = (b-a)*(f(a)+f(b))/2

    # Now evaluate at the other points
    for k = 2:n
        pow2km1 = 2^(k-1);
        h = (b-a)/pow2km1;
        rtable[k,1] = rtable[k-1,1]/2 + h*sum(f(a+i*h) for i in 1:2:pow2km1);
    end

    # Now construct the Romberg table
    for k = 2:n
        for i = k:n
            rtable[i,k] = (4^(k-1)*rtable[i,k-1]-rtable[i-1,k-1])/(4^(k-1)-1);
        end
    end
    rtable
end


end # module
