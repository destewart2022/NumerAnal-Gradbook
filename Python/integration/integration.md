# Integration methods in Python

We have the following integration methods in [`integration.py`](integration.py) in one dimension:

* rectangle rule `rectangle`
* trapezoidal rule `trapezoidal`
* Simpson's rule `simpson`
* mid-point rule `midpt`

and a general composite rule `gen_rule` that uses the weights and nodes for a rule on the interval [0,1].

There is a notebook [`python-integration.ipynb`](python-integration.ipynb) that shows these methods in operation.
