In Kohn-Sham (KS) theory, the ground-state energy of an electronic system (in [[Atomic Units]]) is $$E[\rho] = T_{\text{s}}[\rho] + V_{{\text{ext}}} + J[\rho] + E_{\text{xc}}[\rho],$$where $E_{\text{xc}}[\rho]$ is the [[Exchange-Correlation Energy]], $V_{{\text{ext}}}$ is the [[Electron-Nuclear Energy]], $J[\rho]$ is the [[Hartree Energy]], and $T_{\text{s}}[\rho]$ is the [[Kinetic Energy]]of noninteracting electrons with [[Electron Density]]$$\rho(\textbf{r}) = \sum_{i=1}^N |\varphi_i(\textbf{r})|^2.$$The KS orbitals $\{\varphi_i(\textbf{r})\}$ are solutions of KS equations $$\left\{-\frac{1}{2}\nabla^2 + v_{\text{KS}}([\rho];\textbf{r})\right\}\; \varphi_i(\textbf{r}) = \varepsilon_i\varphi_i(\textbf{r}).$$


## Numerical Quadrature Techniques to Handle the Exchange-Correlation Potential
When the Kohn–Sham matrix is constructed, the exchange-correlation matrix elements $$\begin{equation}V^{\text{xc}}_{\mu\nu}=\int\phi_{\mu}(\textbf{r})V_{\text{xc}}(\textbf{r})\phi_{\nu}(\textbf{r})d\textbf{r}\end{equation}.$$However, in contrast with the other parts of the KS matrix (e.g., Core, electron-electron repulsion), analytical expressions for $\mathbf{V}_{\text{xc}}$ are generally not available. Therefore, [[Numerical Integration|numerical integration]] techniques are to be employed [1]. 

The numerical integration approximates elements of $\mathbf{V}_{\text{xc}}$ as a sum of $P$ points $$\begin{equation}V^{\text{xc}}_{\mu\nu}\approx\hat{V}^{\text{xc}}_{\mu\nu}=\sum^P_{p=1}\phi_{\mu}(\textbf{r}_p)V_{\text{xc}}(\textbf{r}_p)\phi_{\nu}(\textbf{r}_p)\omega_p\end{equation}$$at each point $\textbf{r}_p$ on the grid, weighted by $\omega_p$ whose value depends on the actual numerical technique used [1].
### Implementation
In current quantum chemistry programs, the implementation of the numerical evaluation of the exchange-correlation potential is broken into separate overlapping atomic contributions. To simplify the notation, let us say that $I$ is a value of the integral over an integrand $F(\textbf{r})$$$I = \int F(\textbf{r})d\textbf{r},$$is decomposed into atomic contributions, $A$, over the $M$ nuclei$$I = \sum_A^MI_A,\quad I_A=\int F_A(\textbf{r})d\textbf{r}.$$The atomic integrands $F_A$ are chosen such that heir sum over all nuclei returns the original integrand$$\sum_A^MF_A(\textbf{r})=F(\textbf{r}),$$where the individual $F_A(\textbf{r})$ are constructed from the original integrand by the introduction of
weight functions $\omega_A(\textbf{r})$ with which $F(\textbf{r})$ is multiplied$$F_A(\textbf{r})=\omega_A(\textbf{r})F(\textbf{r}),$$where $\omega_A(\textbf{r})$ assumes a value close to unity if $\textbf{r}$ is close to nucleus $A$ and close to zero near
all other nuclei $B \neq A$. 




# References
1. W. Koch, M.C. Holthausen, *A Chemist's Guide to Density Functional Theory* (John Wiley & Sons, Ltd, 2001), pp. 93-116.