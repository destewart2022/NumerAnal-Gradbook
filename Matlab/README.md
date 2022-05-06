# Matlab code for *Numerical Analysis: A graduate course*

These are Matlab codes for various tasks, following the pseudo-code in the above textbook.

Note that there may be differences due to
* Matlab's 1-relative indexing
* Use of "tracing" or other debugging aids
* Efficient use of function values
* variations in programming/implementation style

### [Solving equations](solvers)
* bisect.m -- bisection method on an interval
* newton.m -- Newton's method (single variable & multivariable)
* gnewton.m -- guarded Newton method (single variable & multivariable)
* secant.m -- secant (derivative free) method
* regfalsi.m -- Regula Falsi method on an interval
* hybrid1.m -- a simple hybrid method

### Interpolation
* lagrange_interp.m -- computes Lagrange interpolation polynomials
* divdif.m -- computes divided differences suitable for interpolation
* dd_interp.m -- computes polynomial interpolant using divided difference table from divdif.m
* dd_hermite.m -- creates divided difference table for Hermite polynomial interpolation (interpolates 1st derivative values as well as function values)
* dd_hinterp.m -- computes Hermite interpolant using results of dd_hermite.m

### Integration
These methods assume that the integrand can accept vectors of evaluation points and produce the vector of function values for the given evaluation points.
* rectangle.m -- composite rectangle rule for integration over an interval
* trapezoidal.m -- composite trapezoidal rule for integration over an interval
* simpson.m -- composite trapezoidal rule for integration over an interval
* gen_int_rule.m -- general integration method for specified weights and nodes for [0,1] to integrate over a given interval
* adapt_simp.m -- adaptive Simpson's rule to integrate a function over an interval with specified accuracy (also uses adapt_simp1.m)
* romberg.m -- Romberg's method using repeated extrapolation to estimate the integral over an interval

### ODE initial value problem solvers
Note that these are all fixed step size, rather than adaptive step-size methods. Using a fixed step size is not only useful for testing things like the order of accuracy of a method, but also for other tasks such as equation solving and optimization as adaption methods can introduce "noise" into function evaluations if that evaluation depends on the result of the ODE solver. 
* euler.m -- Euler's method (explicit)
* heun.m -- Heun's method
* rk4.m -- the standard 4th order Runge-Kutta method
* ab3.m -- the 3rd order Adams-Bashforth method (explicit)
* am3.m -- the 3rd order Adams-Moulton method (implicit)
* More to add!

### Linear algebra routines
TBA!

### Optimization routines
TBA!

