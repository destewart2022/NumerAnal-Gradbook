# Interpolation & approximation
This directory has some Matlab codes for interpolation and approximation, mainly polynomial interpolation

For polynomial interpolation we have
* Lagrange interpolation polynomial evaluations [`lagrange_interp.m`](lagrange_interp.m)
* divided difference evaluation [`divdif.m`](divdif.m) and used to compute the interpolant [`dd_interp.m`](dd_interp.m)
* Hermite interpolation routines [`dd_hermite.m`](dd_hermite.m) and [`dd_hinterp.m`](dd_hinterp.m)

Quasi-interpolation can be carried out using
* Bernstein polynomials on [0,1], which can be computed using [`bernstein.m`](bernstein.m)
