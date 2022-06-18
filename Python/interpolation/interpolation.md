# Interpolation and approximation methods in Python

The file [`interpolation.py`](interpolation.py) contains the functions `divdif` for computing divided differences and `dd_interp` to compute polynomial interpolant using the divided difference table produced by `divdif`.

There is also the notebook [`python-interpolation.ipynb`](python-interpolation.ipynb) that contains some tests and examples of the use of the functions mentioned above.

Cubic spline interpolation is available in Python using [`scipy.interpolate`](https://docs.scipy.org/doc/scipy/reference/interpolate.html). All the usual versions of cubic spline interpolation (natural, clamped, periodic, and not-a-knot) are all available.
