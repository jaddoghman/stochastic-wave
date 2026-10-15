# Linear Stochastic Wave Equation Solvers 

## Overview

The one-dimensional linear stochastic wave equation driven by additive noise and subject to homogeneous Dirichlet boundary conditions is numerically approximated using finite element or finite difference methods for spatial discretization and the explicit stochastic trigonometric integrator for time discretization [1].
The underlying partial differential equation reads

$$
\begin{cases}
  du = vdt, \\
  dv - \Delta u dt = dW,
\end{cases}
$$


## Reference

[1] D. Cohen, On the numerical discretisation of stochastic oscillators, Math. Comput. Simul., (2012). doi:10.1016/j.matcom.2012.02.004.
