# Julia code & notebooks for *Numerical Analysis: A graduate course*

This directory contains Julia code and Julia (.ipynb) notebooks for the above textbook.

The implementation techniques tend to lean on Julia's parametric type system. For example, the trapezoidal method is
declared as `trapezoidal(f::Function,a::T,b::T,n::Int) where {T}` with `T` representing the type of (usually floating point) numbers used
in the computation. This enables the same function to be used for single percision, double precision, or extremely high precision using `BigFloat`s for example, depending on the inputs to `trapezoidal`.
