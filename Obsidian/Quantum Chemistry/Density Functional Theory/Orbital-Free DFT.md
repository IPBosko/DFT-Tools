Orbital-free [[Density Functional Theory|density functional theory]] (OF-DFT) is a practical implementation of the [[Hohenberg–Kohn Theorems|Hohenberg–Kohn theorems]] (HK). This implementation is an alternative to the most popular [[Kohn–Sham DFT]] (KS-DFT).

In OF-DFT theory, only [[Density Functional Approximation|functionals]] that are explicitly dependent on the total [[Electron Density|density]] are used, making it theoretically capable of scaling to much larger systems than KS. The solution of the OF-DFT problem can also be cast in a form which is easy to implement reusing existing KS codes. Within this approach, one solves a KS-like differential equation for only one orbital describing the full system [1]. 

Because of the difficulties in convergence and implementation, *all-electron* implementations of OF-DFT have only been used to derive the energies of small systems, such as atoms and dimers [2, 3]. An all-electron method, refers to methods that solve the [[Schrödinger Equation|Schrödinger equation]] for an electronic system directly in the presence of the $1/r$ non-modified [[Electron-Nuclear Energy|nuclear potential]] and all-electron values are those values that would be obtained with all-electron methods [1].

## Implementation

## References
1. J. Lehtomäki, I. Makkonen, M. A. Caro, A. Harju, O. Lopez-Acevedo, _J. Chem. Phys._ **2014**, _141_, 234102.
2. R. Parr and W. Yang, *Density-functional Theory of Atoms and Molecules*, International Series of Monographs on Chemistry Vol. 16 (Oxford University Press, USA, 1989).
3. V. Karasiev, S. Trickey, _Comput. Phys. Commun._ **2012**, _183_, 2519.


