import numpy as np
from sys import maxsize
import time
from conv_rate_utils import spatial_convergence_rate, plot_convergence_rate
from types import SimpleNamespace

np.set_printoptions(threshold=maxsize, linewidth=maxsize)


# ==========================================================================================
# Stochastic wave equation finite difference solver
# ==========================================================================================

class StoWaveFDSolver:
    def __init__(self, T: float, L: float, Nt: int, Nx: int, s: float, K: int):
        assert s >= 0., "Noise must be either correlated (s>0) or cylindrical (s=0)"
        assert K >= Nx-2, "Number of Wiener modes must be >= number of non-zero dofs"
        self.T, self.L, self.Nt, self.Nx, self.s, self.K = T, L, Nt, Nx, s, K
        dt = T / float(Nt-1)  # Time step size
        self.x = np.linspace(0, L, Nx, endpoint=True)
        self.dx = L / float(Nx-1)  # Spatial step size
        interior_nodes = np.arange(1, Nx-1)

        # Eigenpairs of the Dirichlet Laplacian -Δ
        lam = (4./self.dx**2) * np.sin(interior_nodes*np.pi*self.dx / (2.*L))**2
        phi = np.sqrt(2.*self.dx/L) * np.sin(np.pi*np.outer(self.x[1:-1], interior_nodes) / L)

        # Sine and cosine matrices
        lam_sqrt = np.sqrt(lam)
        self.cos_mat = (phi @ np.diag(np.cos(dt*lam_sqrt))) @ phi.T
        self.sin_mat = (phi @ np.diag(np.sin(dt*lam_sqrt)/lam_sqrt)) @ phi.T

        # Wiener increment components
        modes = np.arange(1, K+1)  # Number of Wiener modes (truncated sum)
        lam_cont_sqrt = modes*np.pi / L  # Square root of continuous eigenvalues of -Δ
        q_sqrt = lam_cont_sqrt**(-s)  # Square root of the eigenvalues of Q
        phi_cont = np.sqrt(2./L) * np.sin(np.pi*np.outer(self.x[1:-1], modes)/L)  # Continuous eigenbasis of -Δ
        self.noise_mat = np.sqrt(dt) * (phi_cont * q_sqrt[None, :] if s > 0 else phi_cont)  # (Nx-2, K)

        # Initial conditions (excluding null boundary conditions)
        self.u0 = np.sin(np.pi*self.x[1:-1]/L)
        self.v0 = np.zeros_like(self.u0, dtype=float)

    def solve(self, xi: np.ndarray):
        """xi: (K, Nt-1) iid N(0,1) random variable generated once for each Monte Carlo sample"""
        dW = self.noise_mat @ xi[:self.K, :]  # Wiener increment (Nx-2, Nt-1)
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

def compute_l2_error_sq(ref: StoWaveFDSolver, coarse: StoWaveFDSolver, u_ref: np.ndarray, u_coarse: np.ndarray):
    """Compute squared L2 error between reference and coarse solutions for convergence rate analysis"""
    Nx_ref = ref.Nx
    Nx_coarse = coarse.Nx
    assert (Nx_ref-1) % (Nx_coarse-1) == 0, "Meshes must be nested"
    stride = (Nx_ref-1) // (Nx_coarse-1)  # To evaluate coarse and reference solutions on the same grid points
    error = coarse.dx * np.sum((u_ref[::stride] - u_coarse)**2)

    return error


# ==========================================================================================
# Application
# ==========================================================================================

if __name__ == "__main__":
    T = 1.  # Final time
    L = 1.  # Interval length
    Nt = 1024  # Number of time grid points
    Nx_values = np.array([16, 32, 64, 128, 4096], dtype=int) + 1  # Last value for reference solution
    M = 1000  # Monte Carlo samples
    s = 1.  # Noise decay rate

    solver_tools = SimpleNamespace(wave_solver=StoWaveFDSolver, l2_squared_error=compute_l2_error_sq)
    start = time.time()
    dx_values, errors = spatial_convergence_rate(solver_tools, T, L, Nt, Nx_values, M, s, parallel=True)
    elapsed = time.time() - start
    print("Elapsed solving time: ", elapsed)

    order = plot_convergence_rate(dx_values, errors, s, method_name="FD")  # Theoretical order
    rate = np.polyfit(np.log(dx_values), np.log(errors), 1)[0]  # Estimated order
    print("Theoretical order: (sup)", order)
    print("Estimated rate: ", rate)
