import numpy as np
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
from os import cpu_count
from threadpoolctl import threadpool_limits
import matplotlib.pyplot as plt
from types import SimpleNamespace


# ==============================================================================
# Common FEM and FD functions
# ==============================================================================

# Sequential mode
# -----------------------------------------------------------------------------
def _spatial_convergence_rate_seq(wave_tools: SimpleNamespace, T: float, L: float, Nt: int, Nx_values: np.ndarray,
                                  M: int, s: float, seed: int):
    # Extract wave_tools variable
    wave_solver = wave_tools.wave_solver  # Class
    compute_l2_error_sq = wave_tools.l2_squared_error  # Function

    # Initializations
    print("[1/2] Initializing reference and coarse solvers...")
    Nx_ref = Nx_values[-1]  # Reference solution mesh size
    K = Nx_ref-2  # Number of Wiener modes
    rng = np.random.default_rng(seed)
    ref = wave_solver(T, L, Nt, Nx_ref, s, K)  # Reference solver
    coarse = np.array([wave_solver(T, L, Nt, nx, s, K) for nx in Nx_values[:-1]])  # Solvers on different meshes
    dx_values = np.array([L / float(nx-1) for nx in Nx_values[:-1]], dtype=float)
    print("[1/2] Initialization complete.")

    # Monte Carlo simulation
    print("[2/2] Solving...")
    errors = np.zeros(coarse.size, dtype=float)
    for _ in tqdm(range(M), desc="Monte Carlo iterations"):
        xi = rng.normal(size=(K, Nt-1))
        u_ref = ref.solve(xi)  # Reference solution
        for i, solver in enumerate(coarse):
            u_coarse = solver.solve(xi)  # Solution on coarse mesh
            errors[i] += compute_l2_error_sq(ref, solver, u_ref, u_coarse)
    errors /= float(M)
    print("[2/2] Solve complete. \n")

    return dx_values, np.sqrt(errors)

# Parallel mode
# -----------------------------------------------------------------------------
_args = {}

def _init_worker(L: float, Nt: int, Nx_values: np.ndarray, ref, coarse, wave_tools: SimpleNamespace):
    threadpool_limits(1)  # Limit threads number used by NumPy to 1
    _args.update(L=L, Nt=Nt, Nx_values=Nx_values, ref=ref, coarse=coarse, wave_tools=wave_tools)

def _single_path_convergence_rate(child_seed):
    rng = np.random.default_rng(child_seed)
    L, Nt, Nx_values = _args["L"], _args["Nt"], _args["Nx_values"]
    ref, coarse = _args["ref"], _args["coarse"]
    wave_tools = _args["wave_tools"]

    Nx_ref = Nx_values[-1]  # Reference solution mesh size
    K = Nx_ref-2  # Number of Wiener modes
    xi = rng.normal(size=(K, Nt-1))
    u_ref = ref.solve(xi)  # Reference solution

    path_errors = np.zeros(coarse.size, dtype=float)  # Pathwise errors
    for i, solver in enumerate(coarse):
        u_coarse = solver.solve(xi)  # Solution on coarse mesh
        path_errors[i] += wave_tools.l2_squared_error(ref, solver, u_ref, u_coarse)

    return path_errors

def _spatial_convergence_rate_par(wave_tools: SimpleNamespace, T: float, L: float, Nt: int, Nx_values: np.ndarray,
                                  M: int, s: float, seed: int):
    # Define random variable seeds
    seed_seq = np.random.SeedSequence(seed)  # Master seed sequence
    child_seeds = seed_seq.spawn(M)  # M independent child seed sequences

    # Convergence rate components
    print("[1/2] Initializing reference and coarse solvers...")
    Nx_ref = Nx_values[-1]  # Reference solution mesh size
    K = Nx_ref-2  # Number of Wiener modes
    ref = wave_tools.wave_solver(T, L, Nt, Nx_ref, s, K)  # Reference solver
    coarse = np.array([wave_tools.wave_solver(T, L, Nt, nx, s, K) for nx in Nx_values[:-1]])  # Solvers on ≠ meshes
    dx_values = np.array([L / float(nx - 1) for nx in Nx_values[:-1]], dtype=float)
    print("[1/2] Initialization complete.")

    print("[2/2] Solving...")
    errors = np.zeros(coarse.size, dtype=float)
    with ProcessPoolExecutor(initializer=_init_worker, initargs=(L, Nt, Nx_values, ref, coarse, wave_tools)) as executor:
        results = executor.map(_single_path_convergence_rate, child_seeds, chunksize=max(1, M // (4*cpu_count())))
        for path_errors in tqdm(results, total=M, desc="Monte Carlo iterations"):
            errors += path_errors
    errors /= float(M)
    print("[2/2] Solve complete. \n")

    return dx_values, np.sqrt(errors)


# ==============================================================================
# Wrapper
# ==============================================================================

def spatial_convergence_rate(wave_tools: SimpleNamespace, T: float, L: float, Nt: int, Nx_values: np.ndarray, M: int,
                             s: float, seed: int=43, parallel: bool=True):
    args = (wave_tools, T, L, Nt, Nx_values, M, s, seed)

    return _spatial_convergence_rate_par(*args) if parallel else _spatial_convergence_rate_seq(*args)


# ==============================================================================
# Plot convergence rate
# ==============================================================================

def plot_convergence_rate(dx_values: np.ndarray, errors: np.ndarray, s: float, method_name: str):
    dim = 1  # Spatial dimension
    beta_max = min(1.+s-dim/2., 3.)
    order = 2.*beta_max/3.  # Theoretical order
    reference = errors[0] * (dx_values / dx_values[0])**order
    plt.figure(figsize=(8, 5))
    plt.loglog(dx_values, errors, marker="o", label="Error " + method_name)
    plt.loglog(dx_values, reference, linestyle="dashed", label=f"Slope {order:.2f}")
    plt.xlabel(r"$\Delta x$ - logarithmic scale", fontsize=14)
    plt.ylabel(r"$L^2(\Omega; L^2)$ error - logarithmic scale", fontsize=14)
    plt.grid(True, which="both")
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.show()

    return order
