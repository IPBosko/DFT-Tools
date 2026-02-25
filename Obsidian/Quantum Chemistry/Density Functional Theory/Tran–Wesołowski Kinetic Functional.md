The [[Kinetic Energy Functional|kinetic energy functional]] developed by Tran and Wesołowski in 2002 [1].
## Overview
This approximate kinetic energy functional utilizes the [[Generalized Gradient Approximation|generalized gradient approximation]](GGA) and is based on the [[Conjointness Conjecture|conjointness conjecture]] proposed by Lee, Lee, and Parr [2]. 
## Mathematical Formulation
The functional adopts the standard GGA form:$$T^{\text{TW}}_{\text{s}}[\rho]=\frac{3}{10}\left(3\pi^2\right)^{2/3}\int\rho^{5/3}(\textbf{r})F^{\text{TW}}_{\text{t}}(s)d\textbf{r}.$$The [[Enhancement Factor|enhancement factor]] $F^{\text{TW}}_{\text{t}}(s)$, where $s$ is the [[Reduced Density Gradient|scaled density gradient]], shares its analytical form with the [[Becke Exchange Functional|Becke86A]] and [[Perdew–Burke–Ernzerhof Exchange-Correlation Functional|PBE]] exchange functionals:$$F^{\text{TW}}_{\text{t}}(s)=1+\kappa-\frac{\kappa}{1+\frac{\mu}{\kappa}s^2}.$$Corresponding Tran–Wesołowski [[Pauli Energy|Pauli functional]]:$$T_{\text{P}}^{\text{TW}}[\rho] = T^{\text{TW}}_{\text{s}}[\rho] - T_{\text{W}}[\rho],$$where $T_{\text{W}}[\rho]$ is the [[von Weizsäcker Kinetic Functional|von Weizsäcker functional]]. Analogously, through the Pauli enhancement factor:$$F^{\text{TW}}_{\text{P}}(s)=1+\kappa-\frac{\kappa}{1+\frac{\mu}{\kappa}s^2} - \frac{5}{3}s^2.$$
## Parametrization
The functional depends on two free parameters, $\kappa$ and $\mu$. Parameters were determined by reproducing the exact kinetic energies of rare-gas atoms. Different fitting sets were used, hence there are four versions of this functional listed on [[Libxc Library|Libxc library]], for details see the original paper [1]. 
## Performance 
In the original paper [1], the functional was tested to reproduce the kinetic energies of twelve closed-shell atoms in [[Self-Consistent-Field Method|post-SCF]] manner. In these tests it outperforms functionals based on the [[von Weizsäcker Kinetic Functional|Gradient Expansion Approximation]] (GEA) and is identified as one of the most accurate approximations despite its simple analytical form.
## References
1. [Tran.2002.IJQC.89.441]
2. [Lee.1991.PRA.44.768]