Most definite integrals cannot be calculated exactly. Value of an integral is approximated numerically when it can't be computed exactly. Another motivation for the numerical integration techniques is that sometimes an algebraic representation of a function is unavailable, instead one has some values of a function at a collection of points.

Since a definite integral is the limit if a sequence of [[Riemann Sum|Riemann sums]], any Riemann sum serves as an approximation of the integral:$$\int_a^bf(x)dx \approx \sum^n_{i=1}f(c_i)\Delta x,$$where $c_i$ is any point chosen from the subinterval $[x_{i-1},x_i]$, for $i=1,2,\dots,n.$ From the definition of definite integral, observe that the larger $n$ is, the better the approximation tends to be. The reason that Riemann sums provide us with numerous approximation schemes is that we are free to choose the evaluation points, $c_i$, for $i=1,2,\dots,n.$ Next section summarizes some of the possible choices of $c_i$.
## Numerical Integration Techniques
There are multiple numerical integration techniques: [[Midpoint Rule|midpoint rule]], [[Trapezoidal Rule|trapezoidal rule]], ...


# References
1. R.T. Smith, R.B. Minton, *Calculus* (McGraw-Hill Science/Engineering/Math, 2003).