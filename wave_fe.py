import numpy as np
from sys import maxsize
from scipy.linalg import solve_banded
import time
from conv_rate_utils import spatial_convergence_rate, plot_convergence_rate
from types import SimpleNamespace

np.set_printoptions(threshold=maxsize, linewidth=maxsize)


# --------------------------------------------------------------------------------
# P1 finite element on 1D uniform mesh (interior nodes only)
# --------------------------------------------------------------------------------

class FiniteElements:
    def __init__(self, x: np.ndarray, n_quad: int=5):
        assert x.ndim == 1 and x.size >= 3, "x must be a 1D array with at least 3 nodes"
        self.x = x
        self.L = x[-1]
        self.Nx = x.size
        self.h = self.L / float(self.Nx-1)
        assert np.allclose(np.diff(self.x), self.h), "Mesh must be uniform"

        # Gauss-Legendre rule mapped to every element
        t, w = np.polynomial.legendre.leggauss(n_quad)
        self._y = x[:-1, None] + 0.5*self.h*(1.+t[None, :])
        self._w = 0.5*self.h * w  # Change of variable constant
        self.phi_right = 0.5 * (1.-t)  # on [x_i, x_{i+1}] (right element)
        self.phi_left = 0.5 * (1.+t)  # on [x_{i-1}, x_i] (left element)

        # Mass matrix reformulation to be used in solve_banded
        self.M_spread = np.zeros((3, self.Nx-2), dtype=float)
        self.M_spread[0, 1:] = self.h / 6.
        self.M_spread[1, :] = 2. * self.h / 3.
        self.M_spread[2, :-1] = self.h / 6.

    def load_vector(self, f):
        """Fill the vector b in Mf_h = b, with b_i = (f, φ_i)"""
        fy = f(self._y)
        if fy.shape != self._y.shape:
            fy = np.vectorize(f, otypes=[float])(self._y)
        b_left = (self._w * fy[:-1, :] * self.phi_left).sum(axis=1)
        b_right = (self._w * fy[1:, :] * self.phi_right).sum(axis=1)

        return b_left + b_right

    def solve_mass_system(self, b: np.ndarray):
        """Solve Mf_h=b for f_h, if b is a known vector"""

        return solve_banded((1, 1), self.M_spread, b)

    def project(self, f):
        """L2 projection of f onto the affine space i.e. solving Mf_h = b for f_h"""
        b = self.load_vector(f)

        return self.solve_mass_system(b)

    def mass_vec_mult(self, vec: np.ndarray):
        """Compute M @ vec in a cheap way"""
        out = 4. * vec
        out[1:] += vec[:-1]
        out[:-1] += vec[1:]

        return (self.h/6.) * out

    def compute_l2_error_sq(self, vec: np.ndarray):
        """Compute L2 error of vec, where vec is a FE vector"""

        return np.dot(vec, self.mass_vec_mult(vec))


# ================================================================================
# Stochastic wave equation finite element solver
# ================================================================================

class StoWaveFEMSolver:
    def __init__(self, T: float, L: float, Nt: int, Nx: int, s: float, K: int):
        assert s >= 0., "Noise must be either correlated (s>0) or cylindrical (s=0)"
        assert K >= Nx-2, "Number of Wiener modes must be >= FE space dimension"
        self.T, self.L, self.Nt, self.Nx, self.s, self.K = T, L, Nt, Nx, s, K
        self.dt = T / float(Nt-1)  # Time step size
        self.x = np.linspace(0, L, Nx, endpoint=True)
        self.fe = FiniteElements(self.x)  # Finite element variable
        h = L / float(Nx-1)  # Spatial step size
        interior_nodes = np.arange(1, Nx-1)

        # Eigenpairs of the Dirichlet Laplacian
        theta = np.pi*h*interior_nodes / L
        lam = (6./h**2) * (1.-np.cos(theta)) / (2.+np.cos(theta))
        phi = np.sqrt(2.*h/L) * np.sin(np.pi*np.outer(self.x[1:-1], interior_nodes)/L)

        # Sine and cosine matrices
        lam_sqrt = np.sqrt(lam)
        self.cos_mat = (phi @ np.diag(np.cos(self.dt*lam_sqrt))) @ phi.T
        self.sin_mat = (phi @ np.diag(np.sin(self.dt*lam_sqrt)/lam_sqrt)) @ phi.T

        # Wiener increments
        modes = np.arange(1, K+1)  # Number of Wiener modes (truncated sum)
        lam_cont_sqrt = modes*np.pi / L  # Square root of continuous eigenvalues of -Δ
        q_sqrt = lam_cont_sqrt**(-s)  # Square root of the eigenvalues of Q
        eta = modes*np.pi*h/(2.*L)
        coef = q_sqrt * (np.sin(eta)/eta)**2 if s > 0 else (np.sin(eta)/eta)**2
        phi_cont = np.sqrt(2./L) * np.sin(np.pi*np.outer(self.x[1:-1], modes)/L)  # Continuous eigenbasis of -Δ
        self.noise_mat = h*np.sqrt(self.dt) * (phi_cont * coef[None, :])  # (Nx-2, K)

        # Initial conditions (L2 projections onto FE space)
        self.u0 = self.fe.project(lambda y: np.sin(np.pi*y/L))
        self.v0 = self.fe.project(lambda y: np.zeros_like(y))

    def solve(self, xi: np.ndarray):
        """xi: (K, Nt-1) iid N(0,1) random variable generated once for each Monte Carlo sample"""
        W_inc_mat = self.noise_mat @ xi[:self.K, :]  # Wiener increments matrix (Nx-2. Nt-1)
        dW = self.fe.solve_mass_system(W_inc_mat)  # Projected Wiener increments
        SdW = self.sin_mat @ dW

        u_prev = self.u0
        u_cur = self.cos_mat @ u_prev + self.sin_mat @ (self.v0 + dW[:, 0])  # Bootstrap first time-step iteration

        # Iterate over all time steps
        for i in range(1, self.Nt-1):
            u_next = 2.*(self.cos_mat @ u_cur) - u_prev + SdW[:, i]
            u_prev, u_cur = u_cur, u_next

        # Append boundary conditions
        u_cur = np.hstack(([0.], u_cur, [0.]))

        return u_cur

def compute_l2_error_sq(ref: StoWaveFEMSolver, coarse: StoWaveFEMSolver, u_ref: np.ndarray, u_coarse: np.ndarray):
    """Compute square L2 error between reference and coarse solutions for convergence rate analysis"""
    u_coarse_fine = np.interp(ref.x, coarse.x, u_coarse)  # Interpolate u_coarse over reference mesh
    error = ref.fe.compute_l2_error_sq(u_ref - u_coarse_fine)

    return error


# ================================================================================
# Application
# ================================================================================

if __name__ == "__main__":
    T = 1.  # Final time
    L = 1.  # Interval length
    Nt = 1024  # Number of time grid points
    Nx_values = np.array([16, 32, 64, 128, 4096], dtype=int) + 1  # Last value for reference solution
    M = 1000  # Monte Carlo samples
    s = 1.  # Noise decay rate

    solver_tools = SimpleNamespace(wave_solver=StoWaveFEMSolver, l2_squared_error=compute_l2_error_sq)
    start = time.time()
    dx_values, errors = spatial_convergence_rate(solver_tools, T, L, Nt, Nx_values, M, s, parallel=True)
    elapsed = time.time() - start
    print("Elapsed time: ", elapsed)

    order = plot_convergence_rate(dx_values, errors, s, method_name="FEM")  # Theoretical order
    rate = np.polyfit(np.log(dx_values), np.log(errors), 1)[0]  # Estimated order
    print("Theoretical order (sup): ", order)
    print("Estimate rate: ", rate)
