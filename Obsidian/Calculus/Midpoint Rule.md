Is a [[Numerical Integration|numerical integration]] technique: $$\int_a^bf(x)dx \approx \sum^n_{i=1}f(c_i)\Delta x,$$where $c_i$ is the midpoint of the subinterval $[x_{i-1},x_i]$: $$c_i=\frac{1}{2}(x_{i-1}+x_i),\quad\text{for}\;i=1,2,\dots,n.$$To obtain higher accuracy, smaller $\Delta x$ could be used, i.e., approximating using more rectangles.
## Algorithm 
1. Store $f(x)$, $a$, $b$, and $n$. 
2. Compute $\Delta x = \frac{b-a}{n}$.
3. Compute $c_1 = a+\frac{\Delta x}{2}$ and start the sum with $f(c_1)$.
4. Compute the next $c_i = c_{i-1}+\Delta x$ and add $f(c_i)$ to the sum.
5. Repeat step 4 until $i=n$, i.e., perform step 4 a total of $(n-1)$ times.
6. Multiply the sum by $\Delta x$.

# References
1. R.T. Smith, R.B. Minton, *Calculus* (McGraw-Hill Science/Engineering/Math, 2003).