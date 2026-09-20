import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt
from time import time
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from tqdm import tqdm 


'''-------------------- VISUALIZATION FUNCTION AND SAVE --------------------'''

def visualize(full_coords, k_val, MAIN_BEAM_RADIUS, MIN_DIST):
    """Generates the array layout plot on the left and the high-resolution SLL radiation pattern map on the right."""
    
    # Calculate adaptive grid size based on aperture frequency bounds
    grid_size = int(k_val * 32)
    print(f"[SPECTRAL CONTROL] Frequency nu={k_val:.4f}. Computation grid dimensions: {grid_size}x{grid_size}")
    
    # Generate 2D space variables (u, v)
    u_bg = np.linspace(-1.0, 1.0, grid_size)
    v_bg = np.linspace(-1.0, 1.0, grid_size)
    U, V = np.meshgrid(u_bg, v_bg)
    Z = np.zeros_like(U)

    # Execute high-fidelity pattern footprint scan
    print("\nLaunching precision radiation pattern scan (CPU)...")
    for i in tqdm(range(U.shape[0]), desc="Computing 2D SLL topology", unit="row"):
        for j in range(U.shape[1]):
            # Constrain grid evaluation strictly inside the visible region (Klein disk R=1)
            if U[i, j] ** 2 + V[i, j] ** 2 <= 1.0:
                Z[i, j] = compute_intensity(
                    np.array([U[i, j], V[i, j]]), full_coords, k_val
                )
            else:
                Z[i, j] = np.nan

    # Convert linear field intensity
    Z_db = 10 * np.log10(Z + 1e-15)

    # Mask the main beam zone
    Z_copy = Z_db.copy()
    main_beam_mask = (U**2 + V**2) < MAIN_BEAM_RADIUS**2
    Z_copy[main_beam_mask] = -np.inf

    # Locate 2D coordinate index and magnitudes of the highest sidelobe peak
    max_idx = np.nanargmax(Z_copy)
    idx_u, idx_v = np.unravel_index(max_idx, Z_copy.shape)

    indep_u = U[idx_u, idx_v]
    indep_v = V[idx_u, idx_v]
    peak_db = Z_db[idx_u, idx_v]
    legend_text = f"Max SLL: {peak_db:.1f} dB"

    # Initialize two-panel layout subplots
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(16, 7), gridspec_kw={"width_ratios": [1, 1.2]}
    )

    # Left Panel: Draw physical minimum element spacing protective zones
    for i in range(len(full_coords)):
        x_c = full_coords[i, 0]
        y_c = full_coords[i, 1]
        safe_circle = Circle(
            (x_c, y_c),
            (MIN_DIST / k_val) / 2,
            color="gray",
            fill=True,
            alpha=0.15,
        )
        ax1.add_patch(safe_circle)
    
    # Plot normalized element coordinate center anchors
    point_width = min(max(1, (MIN_DIST / k_val) / 10), 10)
    ax1.scatter(full_coords[:, 0], full_coords[:, 1], c="black", s=point_width, zorder=3)
    ax1.add_patch(Circle((0, 0), 1.0, color="gray", fill=False, linewidth=1.5))
    ax1.set_xlim(-1.1, 1.1)
    ax1.set_ylim(-1.1, 1.1)
    ax1.set_aspect("equal")
    ax1.set_title("Physical Radiator Layout ($x_k, y_k$)")
    ax1.set_xlabel("$x / R$")
    ax1.set_ylabel("$y / R$")
    ax1.grid(alpha=0.3)

    # Right Panel: Plot contour map of the spectral SLL space
    ax2.set_xlim(-1.05, 1.05)
    ax2.set_ylim(-1.05, 1.05)
    
    cp = ax2.contourf(U, V, Z_db, levels=100, cmap="inferno", vmin=-50, vmax=0)
    fig.colorbar(cp, ax=ax2, label="Relative Intensity (dB)")

    # Render peak marker
    ax2.plot(indep_u, indep_v, "x", color="white", markeredgewidth=1.0, label=legend_text)

    ax2.set_title(f"Planar Array Pattern (nu={k_val:.2f})")
    ax2.set_xlabel("u")
    ax2.set_ylabel("v")
    ax2.legend(loc="upper right", frameon=True, facecolor="black", labelcolor="white", framealpha=0.8)
    ax2.set_aspect("equal", adjustable="box")
    ax2.grid(alpha=0.2)

    # Draw primary main-beam boundary outline
    main_beam_inset = Circle((0, 0), MAIN_BEAM_RADIUS, color="white", fill=False, linestyle="--", alpha=0.7)
    if k_val < 20: # if k_val > 20 circle will be to small
        ax2.add_patch(main_beam_inset)

    # Save precision figure footprint to disk at 300 DPI
    plt.tight_layout()
    plt.savefig(f'interference_plot_n={len(full_coords)}_r={k_val:.1f}.png', bbox_inches='tight', dpi=300)
    plt.show()


def save_coordinates_to_txt(coords, filename="coordinates_2d.txt"):
    """Saves 2D coordinates to a file in custom JSON-like array format."""
    with open(filename, "w") as f:
        f.write("[\n")
        for i, (x, y) in enumerate(coords):
            # Format high-precision string for each coordinate pair
            line = f"  [{x:.16f}, {y:.16f}]"

            # Add separator for all lines except the trailing one
            if i < len(coords) - 1:
                line += ",\n"
            else:
                line += "\n"

            f.write(line)
        f.write("]\n")

    print(f"Coordinates successfully exported to file: {filename}")

'''----------------------- PSLL CALCULATION FUNCTIONS -----------------------'''

@njit
def compute_intensity(uv, full_coords, k_val):
    """Computes array factor normalized intensity at a specific (u, v) point."""
    u, v = uv
    C = 0.0
    S = 0.0
    
    # Pre-compute global constant scaling factor: 2 * pi * nu
    tp_nu = 2.0 * np.pi * k_val
    n_emitters = len(full_coords)

    for i in range(n_emitters):
        x_i = full_coords[i, 0]
        y_i = full_coords[i, 1]
        
        # Calculate localized exponential phase argument
        p_i = tp_nu * (x_i * u + y_i * v)
        
        C += np.cos(p_i)
        S += np.sin(p_i)
           
    return (C**2 + S**2) / (n_emitters**2)


@njit
def find_local_max_newton(m_start, full_coords, k_val):
    """Executes a 2D Newton-Raphson optimization tracking loop to locate a continuous peak point."""
    curr_u = m_start[0]
    curr_v = m_start[1]
    
    # Pre-compute global constant scaling constants
    tp_nu = 2.0 * np.pi * k_val
    tp_nu_sq = tp_nu**2
    
    for _ in range(NEWRON_ITTER):
        C = 0.0
        S = 0.0
        Cu = 0.0
        Su = 0.0
        Cv = 0.0
        Sv = 0.0
        Cuu = 0.0
        Suu = 0.0
        Cvv = 0.0
        Svv = 0.0
        Cuv = 0.0
        Suv = 0.0
        
        # Accumulate quadrature components over all array radiators
        for i in range(len(full_coords)):
            x = full_coords[i, 0]
            y = full_coords[i, 1]
            p = tp_nu * (x * curr_u + y * curr_v)
            
            cos_p = np.cos(p)
            sin_p = np.sin(p)
            
            C += cos_p
            S += sin_p
            
            # Evaluate first partial derivatives of field components
            Cu += -tp_nu * x * sin_p
            Su +=  tp_nu * x * cos_p
            Cv += -tp_nu * y * sin_p
            Sv +=  tp_nu * y * cos_p
            
            # Evaluate second partial derivatives of field components
            cos_p_sq = tp_nu_sq * cos_p
            sin_p_sq = tp_nu_sq * sin_p
            
            Cuu += -x * x * cos_p_sq
            Suu += -x * x * sin_p_sq
            Cvv += -y * y * cos_p_sq
            Svv += -y * y * sin_p_sq
            Cuv += -x * y * cos_p_sq
            Suv += -x * y * sin_p_sq

        # Formulate the analytical array factor intensity gradient vector
        dI_du = 2.0 * (C * Cu + S * Su)
        dI_dv = 2.0 * (C * Cv + S * Sv)
        
        # Evaluate the continuous analytical Hessian matrix coefficients
        H11 = 2.0 * (Cu**2 + Su**2 + C * Cuu + S * Suu)
        H22 = 2.0 * (Cv**2 + Sv**2 + C * Cvv + S * Svv)
        H12 = 2.0 * (Cu * Cv + Su * Sv + C * Cuv + S * Suv)
        
        # Compute Hessian determinant matrix tracking space bounded limits
        det = H11 * H22 - H12**2
        
        # Zero-division safeguard breakout condition
        if abs(det) < 1e-12:
            break
            
        # Execute analytical 2D Newton-Raphson update tracking step
        delta_u = -(H22 * dI_du - H12 * dI_dv) / det
        delta_v = -(H11 * dI_dv - H12 * dI_du) / det
        
        curr_u += delta_u
        curr_v += delta_v

    # Evaluate absolute final intensity at the located terminal coordinates
    final_uv = np.array([curr_u, curr_v])
    final_intensity = compute_intensity(final_uv, full_coords, k_val)
    
    return final_uv, final_intensity


@njit
def find_local_max_on_circle_newton(phi_start, radius, full_coords, k_val):
    """Executes a specialized 1D Newton-Raphson loop constrained strictly along a circular spectrum perimeter."""
    curr_phi = phi_start

    tp_nu = 2.0 * np.pi * k_val
    tp_nu_sq = tp_nu**2
    num_elements = len(full_coords)

    for _ in range(NEWRON_ITTER):
        # Convert angular variable back to spectral coordinates
        u = radius * np.cos(curr_phi)
        v = radius * np.sin(curr_phi)

        C, S = 0.0, 0.0
        Cu, Su, Cv, Sv = 0.0, 0.0, 0.0, 0.0
        Cuu, Suu = 0.0, 0.0

        # Accumulate quadrature fields over all array radiators
        for i in range(num_elements):
            x = full_coords[i, 0]
            y = full_coords[i, 1]
            p = tp_nu * (x * u + y * v)

            cos_p = np.cos(p)
            sin_p = np.sin(p)

            C += cos_p
            S += sin_p

            Cu += -tp_nu * x * sin_p
            Su +=  tp_nu * x * cos_p
            Cv += -tp_nu * y * sin_p
            Sv +=  tp_nu * y * cos_p

            # Calculate only raw uu-components to optimize execution profile
            cos_p_sq = tp_nu_sq * cos_p
            sin_p_sq = tp_nu_sq * sin_p

            Cuu += -x * x * cos_p_sq
            Suu += -x * x * sin_p_sq

        # Formulate core partial derivatives of field intensity
        dI_du = 2.0 * (C * Cu + S * Su)
        dI_dv = 2.0 * (C * Cv + S * Sv)

        # Synthesize remaining Hessian elements via directional vector transformations
        H11 = 2.0 * (Cu**2 + Su**2 + C * Cuu + S * Suu)
        H22 = 2.0 * (Cv**2 + Sv**2 + (v**2 / (u**2 + 1e-20)) * (C * Cuu + S * Suu))
        H12 = 2.0 * (Cu * Cv + Su * Sv + (v / (u + 1e-20)) * (C * Cuu + S * Suu))

        # Evaluate first derivative: dI / dphi (Angular gradient)
        dI_dphi = u * dI_dv - v * dI_du

        # Evaluate second derivative: d^2I / dphi^2 (Angular curvature)
        d2I_dphi2 = H11 * (v**2) - 2.0 * H12 * u * v + H22 * (u**2) - (u * dI_du + v * dI_dv)

        # Execute 1D angular Newton-Raphson correction update step
        curr_phi -= dI_dphi / d2I_dphi2

    # Map terminal tracking state back to cartesian spectrum coordinates
    final_u = radius * np.cos(curr_phi)
    final_v = radius * np.sin(curr_phi)
    final_uv = np.array([final_u, final_v])

    final_intensity = compute_intensity(final_uv, full_coords, k_val)

    return final_uv, final_intensity

'''----------------------- GRADIENT DESCENT FUNCTIONS -----------------------'''

@njit(parallel=True)
def calculate_peak_sll(full_coords, k_val):
    """
    Finds the peak sidelobe level (SLL).
    Combines hexagonal grid and circles (1D and 2D search).
    """

    step = 0.25 / k_val
    main_beam_radius_sq = MAIN_BEAM_RADIUS**2

    step_u = step
    step_v = step * 0.86602540378 

    v_range = np.arange(-1.0, 1.0 + step_v, step_v)
    # Allocation safeguard block for combined 1D/2D search grids
    max_pts = len(v_range) * (int(2.5 / step_u) + 5) + int(40.0 / step) + 500

    starts_arr = np.zeros((max_pts, 2))
    types_arr = np.zeros(max_pts, dtype=np.int32) # 0 - mesh/2D, 1 - inner, 2 - outer
    radii_arr = np.zeros(max_pts)

    count = 0

    # A. Populate baseline hexagonal 2D grid seed points
    for row_idx in range(len(v_range)):
        v = v_range[row_idx]
        u_offset = (step_u / 2.0) if (row_idx % 2 == 1) else 0.0
        u = -1.0 - u_offset
        while u <= 1.0 + step_u:
            if u**2 + v**2 <= 1.0:
                starts_arr[count, 0] = u
                starts_arr[count, 1] = v
                types_arr[count] = 0
                count += 1
            u += step_u

    # B. Populate inner beam perimeter points for 1D angular search
    r_inner = MAIN_BEAM_RADIUS
    d_phi_inner = step / r_inner
    phi_inner_range = np.arange(0.0, 2.0 * np.pi, d_phi_inner)
    for phi in phi_inner_range:
        starts_arr[count, 0] = phi
        types_arr[count] = 1
        radii_arr[count] = r_inner
        count += 1

    # C. Populate outer aperture boundary points for 1D angular search
    r_outer = 1.0
    d_phi_outer = step / r_outer
    phi_outer_range = np.arange(0.0, 2.0 * np.pi, d_phi_outer)
    for phi in phi_outer_range:
        starts_arr[count, 0] = phi
        types_arr[count] = 2
        radii_arr[count] = r_outer
        count += 1

    # D. Inject inner perimeter anchors for independent 2D refinement
    for phi in phi_inner_range:
        starts_arr[count, 0] = r_inner * np.cos(phi)
        starts_arr[count, 1] = r_inner * np.sin(phi)
        types_arr[count] = 0 # Mark as 2D search
        count += 1

    # E. Inject outer perimeter anchors for independent 2D refinement
    for phi in phi_outer_range:
        starts_arr[count, 0] = r_outer * np.cos(phi)
        starts_arr[count, 1] = r_outer * np.sin(phi)
        types_arr[count] = 0 # Mark as 2D search
        count += 1

    num_points = count

    max_vals = np.zeros(num_points)
    coords_u = np.zeros(num_points)
    coords_v = np.zeros(num_points)

    # Execute highly-parallelized combined 1D/2D local tracking loop via Numba threads
    for i in prange(num_points):
        p_type = types_arr[i]

        if p_type == 0:
            m_start = starts_arr[i]
            peak_uv, peak_val = find_local_max_newton(m_start, full_coords, k_val)
        else:
            phi_start = starts_arr[i, 0]
            radius = radii_arr[i]
            peak_uv, peak_val = find_local_max_on_circle_newton(phi_start, radius, full_coords, k_val)
    
        u_p = peak_uv[0]
        v_p = peak_uv[1]
        dist_sq = u_p**2 + v_p**2

        if (dist_sq <= 1.0 and dist_sq >= main_beam_radius_sq) or types_arr[i] != 0:
            max_vals[i] = peak_val
            coords_u[i] = u_p
            coords_v[i] = v_p
        else:
            max_vals[i] = -1.0

    max_idx = np.argmax(max_vals)
    return max_vals[max_idx], np.array([coords_u[max_idx], coords_v[max_idx]])


@njit
def compute_analytical_gradient(pts, k_val, u_peak):
    """
    Computes the intensity gradient with respect to emitter coordinates.
    pts: array of antenna coordinates (N, 2)
    k_val: value of frequency nu
    u_peak: coordinates of the peak (u_m, v_m)
    Output: array of gradients (N, 2)
    """
    u_m, v_m = u_peak[0], u_peak[1]
    n_emitters = pts.shape[0]
    
    tp_nu = 2.0 * np.pi * k_val
    
    # 1. First pass: compute baseline quadratures C and S
    C = 0.0
    S = 0.0
    
    # Array to cache phases and avoid double calculation
    p = np.zeros(n_emitters)
    
    for k in range(n_emitters):
        p[k] = tp_nu * (pts[k, 0] * u_m + pts[k, 1] * v_m)
        C += np.cos(p[k])
        S += np.sin(p[k])
        
    # 2. Second pass: compute gradient components for each antenna element
    grad = np.zeros((n_emitters, 2))
    
    # Pre-compute common dimensional axis scaling coefficients
    coeff_x = 2.0 * tp_nu * u_m
    coeff_y = 2.0 * tp_nu * v_m
    
    for k in range(n_emitters):
        cos_p = np.cos(p[k])
        sin_p = np.sin(p[k])
        
        # Combined analytical tracking common factor for the k-th radiator
        common_factor = S * cos_p - C * sin_p
        
        grad[k, 0] = coeff_x * common_factor # dI / dx_k
        grad[k, 1] = coeff_y * common_factor # dI / dy_k
    
    grad = grad / (n_emitters**2)

    return grad


@njit
def enforce_min_element_spacing_2d(pts, k_val):
    """
    2D projection routine: clamps antenna elements inside the R=1 circle 
    and repels them if the geometric distance falls below d_min.
    """
    n = len(pts)
    d_min = MIN_DIST / k_val
    
    # 1. Enforce boundary conditions for the R=1 circular aperture limit
    for i in range(n):
        r_sq = pts[i, 0]**2 + pts[i, 1]**2
        if r_sq > 1.0:
            r = np.sqrt(r_sq)
            pts[i, 0] /= r
            pts[i, 1] /= r
            
    # 2. Pairwise repulsion processing loop
    for i in range(n):
        for j in range(i + 1, n):
            dx = pts[i, 0] - pts[j, 0]
            dy = pts[i, 1] - pts[j, 1]
            dist = np.sqrt(dx**2 + dy**2)
            
            if dist < d_min and dist > 1e-12:
                overlap = d_min - dist
                # Compute directional split projection vector
                push_x = (dx / dist) * (overlap / 2.0)
                push_y = (dy / dist) * (overlap / 2.0)
                
                pts[i, 0] += push_x
                pts[i, 1] += push_y
                pts[j, 0] -= push_x
                pts[j, 1] -= push_y
                
    # Re-clamp any elements pushed outside the aperture boundary
    for i in range(n):
        r_sq = pts[i, 0]**2 + pts[i, 1]**2
        if r_sq > 1.0:
            r = np.sqrt(r_sq)
            pts[i, 0] /= r
            pts[i, 1] /= r
                
    return pts


@njit
def generate_taylor_circular_coords(N_total):
    """Generates initial coordinates inside R=1 disk based on Taylor density distribution."""
    p = 1.5 # Density decay steepness factor
    final_coords = np.zeros((N_total, 2))
    
    count = 0
    attempts = 0
    max_attempts = N_total * 1000  # Safeguard boundary for large array scales
    
    while count < N_total and attempts < max_attempts:
        # Generate uniform points inside the circle using sqrt radius scaling
        r = np.sqrt(np.random.rand())
        phi = np.random.rand() * 2.0 * np.pi
        
        u = r * np.cos(phi)
        v = r * np.sin(phi)
        
        # Taylor-like window approximation (Max density at r=0, min at r=1)
        accept_prob = (1.0 - r**2)**p
        
        if np.random.rand() < accept_prob:
            final_coords[count, 0] = u
            final_coords[count, 1] = v
            count += 1
        attempts += 1
        
    # Uniform fallback array population if Taylor generator times out
    while count < N_total:
        r = np.sqrt(np.random.rand())
        phi = np.random.rand() * 2.0 * np.pi
        final_coords[count, 0] = r * np.cos(phi)
        final_coords[count, 1] = r * np.sin(phi)
        count += 1
            
    return final_coords


@njit
def optimize_run_2d(n_internal, k_val, iterations, lr):
    """Executes a single optimization pass for the 2D planar array layout."""
    # 1. Generate initial coordinates based on Taylor window inside R=1 disk
    pts = np.zeros((n_internal, 2))
    pts = generate_taylor_circular_coords(n_internal)
    pts = enforce_min_element_spacing_2d(pts, k_val)

    best_pts = pts.copy()
    min_sll = 1e10

    # Initialize logging steps
    print_step = iterations // 100
    if print_step == 0:
        print_step = 1

    for i in range(iterations):
        current_k = k_val

        # Locate global worst sidelobe peak using combined 1D/2D search
        sll_val, u_peak = calculate_peak_sll(pts, current_k)

        if sll_val < min_sll:
            min_sll = sll_val
            best_pts = pts.copy()

        # Print current progression logs
        if i*PRINT_STEP_MULT % print_step == 0:
            percent = (i * 100) // iterations
            sll_db = 10 * np.log10(sll_val)
            print("Progress:", percent, "%, SLL:", sll_db, "dB")

        # Evaluate continuous analytical gradient vector
        grad = compute_analytical_gradient(pts, current_k, u_peak)

        # Normalize gradient magnitudes using global norm
        gnorm = 0.0
        for k in range(n_internal):
            gnorm += grad[k, 0] ** 2 + grad[k, 1] ** 2
        gnorm = np.sqrt(gnorm)

        if gnorm > 1e-12:
            scale = min(max(1.0 / gnorm, MIN_JUMP), 1.0)
            grad = (grad / gnorm) * scale

        # Execute gradient descent step update
        pts = pts - lr * grad

        # Project coordinates back into the circle and repel elements
        pts = enforce_min_element_spacing_2d(pts, k_val)

    return min_sll, best_pts

'''---------------------- CONFIGURATIONS FROM ARTICLE ----------------------'''

#######################
# N_emitters = 100
# nu_val = 4.5
# MIN_DIST = 0.5
# stage1_iters = 100_000
# lr_start = 1e-2
# - 28.6 dB
#######################
# N_emitters = 200
# nu_val = 5.5
# MIN_DIST = 0.5
# stage1_iters = 100_000
# lr_start = 2e-3
# - 33.2 dB
#######################
# N_emitters = 300
# nu_val = 7.0
# MIN_DIST = 0.5
# stage1_iters = 100_000
# lr_start = 1e-3
# - 33.8 dB
#######################
# N_emitters = 500
# nu_val = 9.0
# MIN_DIST = 0.5
# stage1_iters = 100_000
# lr_start = 1e-3
# - 35.3 dB
#######################
# N_emitters = 1000
# nu_val = 12.5
# MIN_DIST = 0.5
# stage1_iters = 100_000
# lr_start = 1e-3
# - 37.0 dB
#######################
#######################
# N_emitters = 600
# nu_val = 60.0
# MIN_DIST = 2.5
# lr_start = 1e-3
# -21.8 dB
#######################
# N_emitters = 2_000
# nu_val = 28.3 * 6.428
# MIN_DIST = 3.21
# PRINT_STEP_MULT = 1
# NEWRON_ITTER = 5
#######################

'''------------------------------- MAIN PART -------------------------------'''

N_emitters = 100
nu_val = 4.5
MIN_DIST = 0.5
PRINT_STEP_MULT = 1
NEWRON_ITTER = 5

dencity = N_emitters/(np.pi*nu_val**2)
print(f"Dencity {dencity:.4f} emitters/lam**2")

MIN_JUMP = 1/nu_val
# main bram radius should be less than 0.817 
MAIN_BEAM_RADIUS = 0.81 / nu_val #0.61 / nu_val

def main():
    print("--- LAUNCHING 2D OPTIMIZATION PASS ---")
    stage1_iters = 100_000
    lr_start = 1e-2
    time_start = time()

    # Launch optimization engine
    print(f"Starting Stage 1 ({stage1_iters} iterations)...")
    final_sll, final_pts = optimize_run_2d(
        n_internal=N_emitters, 
        k_val=nu_val, 
        iterations=stage1_iters, 
        lr=lr_start
    )

    time_stage1 = time() - time_start
    print(f"Completed in {time_stage1:.0f}s. Best SLL = {10 * np.log10(final_sll+1e-15):.2f} dB")
    final_sll_db = 10 * np.log10(final_sll+1e-15)
    
    print("--- OPTIMIZATION RESULTS SUMMARY ---")
    print(f"Final SLL: {final_sll_db:.2f} dB")
        
    # Render and export final layout tracking plots
    print("\nRendering verification plots...")

    filename = f"coordinates_2d_N={N_emitters}_r={nu_val}lam.txt"
    save_coordinates_to_txt(final_pts, filename)

    visualize(final_pts, nu_val, MAIN_BEAM_RADIUS, MIN_DIST)
    

if __name__ == "__main__":
    main()


