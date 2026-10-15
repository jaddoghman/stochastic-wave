# Linear Stochastic Wave Equation Solvers 

## Overview

The one-dimensional linear stochastic wave equation driven by additive noise and subject to homogeneous Dirichlet boundary conditions is numerically approximated on $[0, T] \times [0, L]$ using finite element or finite difference methods for spatial discretization and the explicit stochastic trigonometric integrator for time discretization [1].
The underlying partial differential equation reads

$$
\begin{cases}
  \begin{aligned}
    &\mathrm{d}u = v\mathrm{d}t, \qquad &&\mbox{in } (0, T) \times [0, L], \\
    &dv - \Delta u \mathrm{d}t = \mathrm{d}W, && \mbox{in } (0, T) \times (0, L), \\
    &u(0, \cdot) = u_{0}, \quad v(0, \cdot) = v_{0}, && \mbox{in } [0, L].
  \end{aligned}
\end{cases}
$$

The noise $W$ is a $Q$-Wiener process with covariance operator $Q = (-\Delta)^{-s}$, $s \geq 0$, allowing its Karhunen-Loeve expansion to be expressed in terms of the Dirichlet Laplacian eigenpairs as follows:

$$
W(t, x) = \sum_{k \geq 1}\sqrt{q_{k}}\beta_{k}(t)\phi_{k}(x), \qquad (t, x) \in [0, T] \times [0, L],
$$

where $`\{\phi_{k}(x)\}_{k \geq 1}`$ is an orthonormal basis of $`L^{2}(0, L)`$ defined by $`\phi_{k}(x) = \sqrt{2/L}\sin(k\pi x/L)`$, $`\{\beta_{k}\}_{k \geq 1}`$ is a family of mutually independent identically distributed real-valued Brownian motions, $`q_{k} = \lambda_{k}^{-s}`$, and $`\{\lambda_{k}\}_{k \geq 1}`$ is the set of eigenvalues of $`-\Delta`$ associated with $`\{\phi_{k}\}_{k \geq 1}`$ defined by $`\lambda_{k} = (k\pi/L)^{2}`$.


## Reference

[1] D. Cohen, On the numerical discretisation of stochastic oscillators, Math. Comput. Simul., (2012). doi:10.1016/j.matcom.2012.02.004.

[2] M. Kovács, S. Larsson, F. Saedpanah, Finite element approximation of the linear stochastic wave equation with additive noise, SIAM J. Numer. Anal., (2010), doi:10.1137/090772241.

[3] D. Cohem, S. Larsson, M. Sigg, A trigonometric method for the linear stochastic wave equation, SIAM J. Numer. Anal., (2013). doi:10.1137/12087030X.
