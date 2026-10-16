# Linear Stochastic Wave Equation Solvers 

## Overview

The one-dimensional linear stochastic wave equation driven by additive noise and subject to homogeneous Dirichlet boundary conditions is numerically approximated on $[0, T] \times [0, L]$ using **finite element** or **finite difference** methods for spatial discretization and the explicit stochastic **trigonometric integrator** for time discretization [1].
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

The noise $W$ is a $Q$-**Wiener process**, cylindrical or correlated, with covariance operator $`Q = (-\Delta)^{-s}`$, $s \geq 0$, allowing its Karhunen-Loeve expansion to be expressed in terms of the Dirichlet Laplacian eigenbasis as follows:

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
  & U^{n} = C_{\Delta x}^{\Delta t}U^{n-1} + \Omega_{\Delta x}^{-1}S_{\Delta x}^{\Delta t}(V^{n-1} + \Delta_{n}W), \\
  & V^{n} = -\Omega_{\Delta x}S_{\Delta x}^{\Delta t}U^{n-1} + C_{\Delta x}^{\Delta t}(V^{n-1} + \Delta_{n}W),
\end{cases}
$$

where $`(\Delta_{n}W)_{i} = W(t_{n}, x_{i}) - W(t_{n-1}, x_{i}) = \sum_{k=1}^{K}\sqrt{q_{k}}\Delta_{n}\beta_{k}\phi_{k}(x_{i})`$, with $`\Delta_{n}\beta_{k} \sim \mathcal{N}(0, \Delta t)`$.
Alternatively, the above fully discrete scheme can be turned into the following two-step recurrence:

$$
\begin{cases}
  & U^{n+1} = 2C_{\Delta x}^{\Delta t}U^{n} - U^{n-1} + \Omega_{\Delta x}^{-1}S_{\Delta x}^{\Delta t}\Delta_{n+1}W, \\
  & U^{1} = C_{\Delta x}^{\Delta t}U^{0} + \Omega_{\Delta}^{-1}S_{\Delta x}^{\Delta t}(V^{0} + \Delta_{1}W)
\end{cases}
$$

for all $`n \in \{1, \dotsc, N_{t}-2\}`$.
For more details on the finite difference scheme, the reader may refer to [4].

### Finite element method

For notation consistency, denote $h = \Delta x$, and $\Lambda_{h} = -\Delta$.
The $P_{1}$ elements are adopted, the corresponding finite element space is of dimension $N_{x}-2$ due to the homogeneous Dirichlet boundary conditions, and the associated mass and stiffness matrices are denoted by $M$ and $K$, respectively.
Denote by $`\{\varphi_{i}\}_{1 \leq i \leq N_{x}-2}`$ (not to be confused with $\phi_{i}$) the set of hat functions defined over $[0, L]$.
The discrete eigenpairs of the approximated Dirichlet Laplacian $\Lambda_{h}$ are given by:

$$
\lambda_{k}^{h} = \frac{6}{h^{2}}\frac{1-\cos(\theta_{k})}{2+\cos(\theta_{k})}, \qquad (\phi_{k}^{h})_{i} = \sqrt{\frac{2h}{L}}\sin(\frac{k\pi x_{i}}{L}),
$$

where $\theta_{k} = k\pi h/L$, for all $`i, k \in \{1, \dotsc, N_{x}-2\}`$.
The coordinate representation of the discrete Laplace operator in the nodal basis is $\Lambda_{h} = M^{-1}K$, but does not require an explicit computation since it is diagonalizable and can be replaced by its diagonal matrix.
By jointly applying the finite element method and the explicit stochastic trigonometric time integrator, the fully discrete scheme reads

$$
\begin{cases}
  & U^{n} = C_{h}^{\Delta t}U^{n-1} + \Omega_{h}^{-1}S_{h}^{\Delta t}(V^{n-1} + M^{-1}\Delta_{n}W), \\
  & V^{n} = -\Omega_{h}S_{h}^{\Delta t}U^{n-1} + C_{h}^{\Delta t}(V^{n-1} + M^{-1}\Delta_{n}W),
\end{cases}
$$

where $\Omega_{h} = \Lambda_{h}$, $C_{h}^{\Delta t} = \cos(\Delta t\Omega_{h})$, $S_{h}^{\Delta t} = \sin(\Delta t\Omega_{h})$, and $`(\Delta_{n}W)_{i} = \sum_{k=1}^{K}\sqrt{q_{k}}\Delta_{n}\beta_{k}\int_{0}^{L}\phi_{k}\varphi_{i}\mathrm{d}x`$, with $\Delta_{n}\beta_{k} \sim \mathcal{N}(0, \Delta t)$.
The discrete noise term $M^{-1}\Delta_{n}W$ is simply the $L^{2}$ projection of $\Delta_{n}W$ onto the finite element space, and does not require forming $M^{-1}$.
The above fully discrete scheme can be recast as a two-step recurrence:

$$
\begin{cases}
  & U^{n+1} = 2C_{h}^{\Delta t}U^{n} - U^{n-1} + \Omega_{h}^{-1}S_{h}^{\Delta t}M^{-1}\Delta_{n+1}W, \\
  & U^{1} = C_{h}^{\Delta t}U^{0} + \Omega_{h}^{-1}S_{h}^{\Delta t}(V^{0} + M^{-1}\Delta_{1}W),
\end{cases}
$$

for all $`n \in \{1, \dotsc, N_{t}-2\}`$.
The reader may refer to [3] for further details.

## Reference

[1] D. Cohen. On the numerical discretisation of stochastic oscillators, Math. Comput. Simul., (2012). doi:10.1016/j.matcom.2012.02.004.

[2] M. Kovács, S. Larsson, F. Saedpanah. Finite element approximation of the linear stochastic wave equation with additive noise, SIAM J. Numer. Anal., (2010), doi:10.1137/090772241.

[3] D. Cohem, S. Larsson, M. Sigg. A trigonometric method for the linear stochastic wave equation, SIAM J. Numer. Anal., (2013). doi:10.1137/12087030X.

[4] D. Cohen, L. Quer-Sardanyons. A fully discrete approximation of the one-dimensional stochastic wave equation, IMAJNA, (2015), doi:10.1093/imanum/drn000
