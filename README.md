# Linear Stochastic Wave Equation Solvers 

## Overview

The one-dimensional linear stochastic wave equation driven by additive noise and subject to homogeneous Dirichlet boundary conditions is numerically approximated on $[0, T] \times [0, L]$ using **finite element** or **finite difference** methods for spatial discretization and the explicit stochastic **trigonometric integrator** for time discretization [1].
The underlying partial differential equation reads

$$
\begin{cases}
  \begin{aligned}
    \mathrm{d}u = v\mathrm{d}t, \qquad &\mbox{in } (0, T) \times [0, L], \\
    dv - \Delta u \mathrm{d}t = \mathrm{d}W, & \mbox{in } (0, T) \times (0, L), \\
    u(0, \cdot) = u_{0}, \quad v(0, \cdot) = v_{0}, & \mbox{in } [0, L].
  \end{aligned}
\end{cases}
$$

The noise $W$ is a **Wiener process**, cylindrical or correlated, with covariance operator $`Q = (-\Delta)^{-s}`$, $s \geq 0$, allowing its Karhunen-Loève expansion to be expressed in terms of the Dirichlet Laplacian eigenbasis as follows:

$$
W(t, x) = \sum_{k \geq 1}\sqrt{q_{k}}\beta_{k}(t)\phi_{k}(x), \qquad (t, x) \in [0, T] \times [0, L],
$$

where $`\{\phi_{k}\}_{k \geq 1}`$ is an orthonormal basis of $`L^{2}(0, L)`$ defined by $`\phi_{k}(x) = \sqrt{2/L}\sin(k\pi x/L)`$, $`\{\beta_{k}\}_{k \geq 1}`$ is a family of mutually independent and identically distributed real-valued Brownian motions, $`q_{k} = \lambda_{k}^{-s}`$, and $`\{\lambda_{k}\}_{k \geq 1}`$ is the set of eigenvalues of $`-\Delta`$ associated with $`\{\phi_{k}\}_{k \geq 1}`$ defined by $`\lambda_{k} = (k\pi/L)^{2}`$.

## Discretization

The time inteval $[0, T]$ is discretized into $`N_{t} \in \mathbb{N}\setminus\{0, 1\}`$ equispaced nodes with time step size $\Delta t = T/(N_{t}-1)$ and grid points $t_{n} = n\Delta t$, for all $`n \in \{0, 1, \dotsc, N_{t}-1\}`$.
The physical domain $[0, L]$ is discretized using $`N_{x} \in \mathbb{N}\setminus\{0, 1\}`$ grid points, with a uniform spatial step size $\Delta x = L/(N_{x}-1)$, and nodes $x_{i} = i\Delta x$, for all $`i \in \{0, 1, \dotsc, N_{x}-1\}`$.

### Finite difference method

The Laplace operator $\Delta$ is discretized using the second-order central difference scheme, resulting in the following discrete eigenpairs of $\Lambda_{\Delta x} = -\Delta$:

$$
\lambda_{k}^{\Delta x} = \frac{4}{\Delta x^{2}}\sin^{2}(\frac{k\pi\Delta x}{2L}), \qquad (\phi_{k}^{\Delta x})_{i} = \sqrt{\frac{2\Delta x}{L}}\sin(\frac{k\pi x_{i}}{L}).
$$

Setting $\Omega_{\Delta x} = \Lambda_{\Delta x}^{1/2}$, $C_{\Delta x}^{\Delta t} = \cos(\Delta t\Omega_{\Delta x})$, and $S_{\Delta x}^{\Delta t} = \sin(\Delta t\Omega_{\Delta x})$, the fully discrete scheme, combining the finite difference method and the explicit stochastic trigonometric time integrator, is given by:

$$
\begin{cases}
  U^{n} = C_{\Delta x}^{\Delta t}U^{n-1} + \Omega_{\Delta x}^{-1}S_{\Delta x}^{\Delta t}(V^{n-1} + \Delta_{n}W), \\
  V^{n} = -\Omega_{\Delta x}S_{\Delta x}^{\Delta t}U^{n-1} + C_{\Delta x}^{\Delta t}(V^{n-1} + \Delta_{n}W),
\end{cases}
$$

where $`(\Delta_{n}W)_{i} = W(t_{n}, x_{i}) - W(t_{n-1}, x_{i}) \approx \sum_{k=1}^{K}\sqrt{q_{k}}\Delta_{n}\beta_{k}\phi_{k}(x_{i})`$, with $K \leq N_{x}-2$ the number of modes and $`\Delta_{n}\beta_{k} \sim \mathcal{N}(0, \Delta t)`$.
Alternatively, the above fully discrete scheme can be turned into the following two-step recurrence:

$$
\begin{cases}
  U^{n+1} = 2C_{\Delta x}^{\Delta t}U^{n} - U^{n-1} + \Omega_{\Delta x}^{-1}S_{\Delta x}^{\Delta t}\Delta_{n+1}W, \\
  U^{1} = C_{\Delta x}^{\Delta t}U^{0} + \Omega_{\Delta x}^{-1}S_{\Delta x}^{\Delta t}(V^{0} + \Delta_{1}W)
\end{cases}
$$

for all $`n \in \{1, \dotsc, N_{t}-2\}`$.
For more details on the finite difference scheme, the reader may refer to [4].

### Finite element method

For notation consistency, denote $h = \Delta x$, and $\Lambda_{h} = -\Delta$.
The $P_{1}$ elements are adopted, the corresponding finite element space is of dimension $N_{x}-2$ due to the homogeneous Dirichlet boundary conditions, and the associated mass matrix is denoted by $M$.
Denote by $`\{\varphi_{i}\}_{1 \leq i \leq N_{x}-2}`$ (not to be confused with $\phi_{i}$) the set of hat functions defined over $[0, L]$.
The discrete eigenpairs of the approximated Dirichlet Laplacian $\Lambda_{h}$ are given by:

$$
\lambda_{k}^{h} = \frac{6}{h^{2}}\frac{1-\cos(\theta_{k})}{2+\cos(\theta_{k})}, \qquad (\phi_{k}^{h})_{i} = \sqrt{\frac{2h}{L}}\sin(\frac{k\pi x_{i}}{L}),
$$

where $\theta_{k} = k\pi h/L$, for all $`i, k \in \{1, \dotsc, N_{x}-2\}`$.
By jointly applying the finite element method and the explicit stochastic trigonometric time integrator, the fully discrete scheme reads

$$
\begin{cases}
  U^{n} = C_{h}^{\Delta t}U^{n-1} + \Omega_{h}^{-1}S_{h}^{\Delta t}(V^{n-1} + M^{-1}\Delta_{n}W), \\
  V^{n} = -\Omega_{h}S_{h}^{\Delta t}U^{n-1} + C_{h}^{\Delta t}(V^{n-1} + M^{-1}\Delta_{n}W),
\end{cases}
$$

where $\Omega_{h} = \Lambda_{h}^{1/2}$, $C_{h}^{\Delta t} = \cos(\Delta t\Omega_{h})$, $S_{h}^{\Delta t} = \sin(\Delta t\Omega_{h})$, and $`(\Delta_{n}W)_{i} \approx \sum_{k=1}^{K}\sqrt{q_{k}}\Delta_{n}\beta_{k}\int_{0}^{L}\phi_{k}\varphi_{i}\mathrm{d}x`$, with $K \geq N_{x}-2$ (see [2]) and $\Delta_{n}\beta_{k} \sim \mathcal{N}(0, \Delta t)$.
The discrete noise term $M^{-1}\Delta_{n}W$ is simply the $L^{2}$ projection of $\Delta_{n}W$ onto the finite element space, and does not require forming $M^{-1}$.
The above fully discrete scheme can be recast as a two-step recurrence:

$$
\begin{cases}
  U^{n+1} = 2C_{h}^{\Delta t}U^{n} - U^{n-1} + \Omega_{h}^{-1}S_{h}^{\Delta t}M^{-1}\Delta_{n+1}W, \\
  U^{1} = C_{h}^{\Delta t}U^{0} + \Omega_{h}^{-1}S_{h}^{\Delta t}(V^{0} + M^{-1}\Delta_{1}W),
\end{cases}
$$

for all $`n \in \{1, \dotsc, N_{t}-2\}`$.
The reader may refer to [3] for further details.

## Solvers

The repository contains finite element (`wave_fe.py`) and finite difference (`wave_fd.py`) solvers encapsulated in `StoWaveFEMSolver` and `StoWaveFDSolver` as Python classes.
To keep this repository self-contained, the $P_{1}$ elements along with some associated features are implemented in `FiniteElements` inside `wave_fe.py`.
The finite element solver can be invoked as follows:

```python
solver_fe = StoWaveFEMSolver(T, L, Nt, Nx, s, K)  # Create fem solver
u_fe = solver_fe.solve(xi)  # Compute fem solution
```

where `K` is the number of (truncated) Wiener modes, and `xi` $\sim \mathcal{N}(0, 1)$ is an array of size $(K, N_{t}-1)$.
The computed finite element solution `u_fe` is a single sample path (realization) and can vary from one simulation to another.

The finite difference solver preserves the same syntax and outputs:

```python
solver_fd = StoWaveFDSolver(T, L, Nt, Nx, s, K)  # Create fd solver
u_fd = solver_fd.solve(xi)  # Compute fd solution
```

## Convergence rate analysis

The spatial convergence rate analysis involves Monte Carlo simulations, which can make the overall computational cost high.
The module `conv_rate_utils.py` implements convergence rate functions in sequential and parallel modes; each function includes both underlying numerical methods to avoid code duplication.
The sequential and parallel implementations are both wrapped in one function:

```python
spatial_convergence_rate(solver_tools, T, L, Nt, Nx_values, M, s, seed, parallel)
```

where `solver_tools` is a SimpleNamespace with attributes `wave_solver` and `l2_squared_error`, `M` is the number of Monte Carlo samples (not to be confused with the mass matrix), `seed` is an integer for reproducible random results, and `parallel` is a boolean to switch between sequential and parallel modes.
This function returns two arrays containing the spatial step size and the $L^{2}(\Omega; L^{2}(0, L))$ error at time $T$, where $\Omega$ is the probability space.
Since the method of manufactured solutions cannot be used to obtain an exact solution, the spatial convergence rate is assessed by comparing a reference solution computed on a fine mesh with other solutions obtained on coarser meshes.

This same module provides a function to plot the computed spatial convergence rate of the displacement solution $u$:

```python
plot_convergence_rate(dx_values, errors, s, method_name)
```

where `dx_values` and `errors` are produced by `spatial_convergence_rate`, `s` is the decay rate of $Q = \Lambda^{-s}$, and `method_name` is a string (e.g., FEM or FD).
In addition to plotting, this function returns the theoretical order, which can be compared with the numerical order using NumPy's `polyfit` function.

### Theoretical results

Given that the stochastic wave equation is implemented in 1D, the finite element and finite difference solvers have equivalent spatial convergence rates (same eigenvectors, both with eigenvalue error $`\mathcal{O}(k^{4}h^{2})`$).
Based on [2], the theoretical $L^{2}$ spatial convergence rate (upper bound) is $2\beta/3$, where $\beta \in [0, 3]$ is the regularity of the solution.
In short, for $s \geq 0$, $\beta < 1+s-d/2$, where $d$ is the spatial dimension.

## Features

* One-dimensional linear stochastic wave equation solvers
* Finite element spatial discretization
* Finite difference spatial discretization
* Explicit stochastic trigonometric time integrator
* Automated computation of $L^{2}$ errors
* Numerical convergence analysis
* Sequential and Parallelized functions
* Validation of theoretical results

## Requirements

* Python 3.x
* NumPy
* Scipy
* threadpoolctl
* tqdm
* matplotlib

Additional packages may be required depending on the plotting and post-processing tools used.

## Running

Run a finite element simulation with

```bash
python3 wave_fe.py
```

or a finite difference simulation with

```bash
python3 wave_fd.py
```

## Purpose

This repository was developed as a numerical study of explicit trigonometric time integration for finite element and finite difference discretizations of hyperbolic partial differential equations.
It is intended for educational and research purposes, particularly for verifying implementations and investigating convergence behavior.

## References

[1] D. Cohen. On the numerical discretisation of stochastic oscillators, Math. Comput. Simul., (2012). doi:10.1016/j.matcom.2012.02.004.

[2] M. Kovács, S. Larsson, F. Saedpanah. Finite element approximation of the linear stochastic wave equation with additive noise, SIAM J. Numer. Anal., (2010), doi:10.1137/090772241.

[3] D. Cohen, S. Larsson, M. Sigg. A trigonometric method for the linear stochastic wave equation, SIAM J. Numer. Anal., (2013). doi:10.1137/12087030X.

[4] D. Cohen, L. Quer-Sardanyons. A fully discrete approximation of the one-dimensional stochastic wave equation, IMAJNA, (2016), doi:10.1093/imanum/drv006
