import numpy as np
import quadpy
from math import pi, sinh, cos, log, sqrt
import h5py
import xml.etree.ElementTree as ET

def get_time_days_from_h5(f, t_idx):
    """
    Read physical time (days) from Mesh_Spatial_Domain_<t_idx> XML stored in the HDF5.
    """
    mesh_key = f"Mesh_Spatial_Domain_{t_idx}"
    if mesh_key not in f:
        raise KeyError(f"Missing dataset '{mesh_key}' in H5 file.")
    root = ET.fromstring(f[mesh_key][:])
    node = root.find("Time")
    return float(node.attrib["Value"])


def _tracy_precompute_params(t_days,
                             a, L, alpha, psi_r,
                             theta_s, theta_r, Ks_day,
                             n_terms):
    # basic parameters
    h0 = 1.0 - np.exp(alpha * psi_r)        # \bar h_0
    Ks = Ks_day / 86400.0                   # m/s
    c  = alpha * (theta_s - theta_r) / Ks   # 1/s
    t  = t_days * 86400.0                   # s
    # eigenvalues and decay rates
    k_idx  = np.arange(1, n_terms + 1, dtype=float)  # k = 1..n_terms
    lam    = k_idx * pi / L                          # λ_k
    gamma1 = (lam**2 + (alpha**2) / 4.0) / c         # γ_1
    gamma2 = ((2.0 * pi / a)**2 + lam**2 + (alpha**2) / 4.0) / c  # γ_2
    # mode coefficients that depend on t but not on x,z
    sign = (-1.0) ** k_idx
    A = sign * lam * (1.0 / gamma1) * np.exp(-gamma1 * t)  # for first part
    B = sign * lam * (1.0 / gamma2) * np.exp(-gamma2 * t)  # for second part
    # steady-state root
    root_ss = sqrt((alpha / 2.0)**2 + (2.0 * pi / a)**2)
    return {
        "h0": h0,
        "Ks": Ks,
        "c": c,
        "t": t,
        "lam": lam,
        "A": A,
        "B": B,
        "root_ss": root_ss,
        "alpha": alpha,
        "psi_r": psi_r,
        "a": a,
        "L": L,
    }
def _tracy_eval(nodes, params):
    """
    Vectorized evaluation of Tracy transient solution at a set of points.

    nodes : (N,2) or (N,3), uses columns 0:(x) and 1:(z)
    params: precomputed dictionary from _tracy_precompute_params
    """
    x = nodes[:, 0]
    z = nodes[:, 1]

    h0      = params["h0"]
    c       = params["c"]
    lam     = params["lam"]
    A       = params["A"]
    B       = params["B"]
    root_ss = params["root_ss"]
    alpha   = params["alpha"]
    psi_r   = params["psi_r"]
    a       = params["a"]
    L       = params["L"]

    psi_tr = np.empty_like(x)

    # mask for top boundary z = L (Dirichlet condition)
    top_mask = np.abs(z - L) < 1e-12
    interior_mask = ~top_mask

    # --- top boundary: direct BC formula ---
    if np.any(top_mask):
        x_top = x[top_mask]
        psi_tr[top_mask] = (1.0 / alpha) * np.log(
            np.exp(alpha * psi_r)
            + 0.5 * h0 * (1.0 - np.cos(2.0 * np.pi * x_top / a))
        )

    # --- interior points: steady + transient ---
    if np.any(interior_mask):
        xi = x[interior_mask]
        zi = z[interior_mask]

        # steady transformed part \bar h_ss(x,z)
        exp_factor = np.exp(0.5 * alpha * (L - zi))
        sinh_alpha_zi = np.sinh(0.5 * alpha * zi)
        sinh_alpha_L  = np.sinh(0.5 * alpha * L)
        cos_x = np.cos(2.0 * np.pi * xi / a)
        sinh_root_zi = np.sinh(root_ss * zi)
        sinh_root_L  = np.sinh(root_ss * L)

        h_bar_ss = (
            0.5 * h0
            * exp_factor
            * (
                sinh_alpha_zi / sinh_alpha_L
                - cos_x * (sinh_root_zi / sinh_root_L)
            )
        )
        # transient correction \bar\phi(x,z,t)
        sin_mat = np.sin(lam[:, None] * zi[None, :])   # (K,Ni)
        SA = (A[:, None] * sin_mat).sum(axis=0)        # (Ni,)
        SB = (B[:, None] * sin_mat).sum(axis=0)        # (Ni,)
        s = SA - cos_x * SB
        phi_bar = (h0 / (L * c)) * exp_factor * s
        # total transformed variable and back-transform
        h_bar = h_bar_ss + phi_bar
        psi_tr[interior_mask] = (1.0 / alpha) * np.log(
            np.exp(alpha * psi_r) + h_bar
        )

    return psi_tr

def tracy_transient(nodes, t_days, *,
                    a, L, alpha, psi_r,
                    theta_s, theta_r, Ks_day,
                    n_terms=120):
    """
    Transient Tracy solution psi(x,z,t) for each node.
    """
    params = _tracy_precompute_params(
        t_days,
        a, L, alpha, psi_r,
        theta_s, theta_r, Ks_day,
        n_terms
    )
    return _tracy_eval(nodes, params)


# ----------------------------------------------------------------------
# L2 / Linf error using quadpy with precomputed Tracy params
# ----------------------------------------------------------------------
def tracy_transient_L2_Linf_error_quadpy_psi(u_num,
                                             nodes,
                                             elements,
                                             *,
                                             t_days,
                                             a, L, alpha, psi_r,
                                             theta_s, theta_r, Ks_day,
                                             n_terms=120,
                                             degree=5):
    """
    Transient Tracy error in psi only: returns (L2_psi, Linf_psi)
    """
    scheme = quadpy.t2.get_good_scheme(degree)

    tracy_params = _tracy_precompute_params(
        t_days,
        a, L, alpha, psi_r,
        theta_s, theta_r, Ks_day,
        n_terms
    )

    L2_total   = 0.0
    Linf_total = 0.0

    for elem in elements:
        coords  = nodes[elem]
        coords2 = coords[:, :2]       # (x,z)
        u_loc   = u_num[elem]         # nodal psi on this element

        v0 = coords2[0]
        J  = np.column_stack((coords2[1] - v0, coords2[2] - v0))
        detJ = np.linalg.det(J)
        if abs(detJ) < 1e-30:
            continue
        invJ = np.linalg.inv(J)

        # vertex Linf
        psi_ex_v = _tracy_eval(coords2, tracy_params)
        local_Linf = float(np.max(np.abs(u_loc - psi_ex_v)))

        def integrand(xq):
            nonlocal local_Linf

            pts = xq.T
            dx  = pts - v0[None, :]
            lam12 = dx @ invJ.T
            lam1  = lam12[:, 0]
            lam2  = lam12[:, 1]
            lam0  = 1.0 - lam1 - lam2

            psi_num_q = lam0*u_loc[0] + lam1*u_loc[1] + lam2*u_loc[2]
            psi_ex_q  = _tracy_eval(pts, tracy_params)

            diff = psi_num_q - psi_ex_q
            local_Linf = max(local_Linf, float(np.max(np.abs(diff))))
            return diff**2

        L2_total += scheme.integrate(integrand, coords2)
        Linf_total = max(Linf_total, local_Linf)

    return np.sqrt(L2_total), Linf_total



# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

Lx = 10.0  # domain length

ref_info = [
#    ("ref_0", 11),
#    ("ref_1", 21),
    ("ref_2", 41),
    ("ref_3", 81),
    ("ref_4", 161),
    ("ref_5", 321),
    ("ref_6", 641)
]

time_indices = [10, 20, 50, 100, 200, 400] #[10, 20, 50, 100, 200]  

OUT_TXT = "Tracy_transient_L2_Linf_convergence_new_4.txt"

with open(OUT_TXT, "w") as fout:

    for t_idx in time_indices:
        h5_path0 = f"{ref_info[0][0]}/re_vgm_sand_10x10m_2d.h5"
        with h5py.File(h5_path0, "r") as f0:
            t_days = get_time_days_from_h5(f0, t_idx)

        fout.write(f"\n===== Errors at t = {t_days:.6e} days (index={t_idx}) =====\n")
        results = []
        for ref_name, nnx in ref_info:
            h5_path = f"{ref_name}/re_vgm_sand_10x10m_2d.h5"
            with h5py.File(h5_path, "r") as f:
                elements = f["elementsSpatial_Domain0"][:]
                nodes    = f["nodesSpatial_Domain0"][:]
                psi_num  = f[f"pressure_head_t{t_idx}"][:]
                #t_days   = get_time_days_from_h5(f, t_idx)
            
            if elements.min() == 1:
                elements = elements - 1

            L2_psi, Linf_psi = tracy_transient_L2_Linf_error_quadpy_psi(
                psi_num,
                nodes,
                elements,
                t_days=t_days,
                a=10.0,
                L=10.0,
                alpha=0.164,
                psi_r=-15.24,
                theta_s=0.301,
                theta_r=0.093,
                Ks_day=2.04,
                n_terms=200,
                degree=5,
            )

            h = Lx / (nnx - 1)
            results.append([ref_name, nnx, h, L2_psi, Linf_psi])

        rates_psi_L2   = [None]
        rates_psi_Linf = [None]

        for k in range(1, len(results)):
            hC, hF = results[k-1][2], results[k][2]

            EC_L2, EF_L2 = results[k-1][3], results[k][3]
            EC_Li, EF_Li = results[k-1][4], results[k][4]

            pL2 = np.log(EC_L2 / EF_L2) / np.log(hC / hF) if EF_L2 > 0 else None
            pLi = np.log(EC_Li / EF_Li) / np.log(hC / hF) if EF_Li > 0 else None

            rates_psi_L2.append(pL2)
            rates_psi_Linf.append(pLi)

        fout.write("\n--- Convergence in psi (pressure head) ---\n")
        fout.write(f"{'ref':<6} {'nnx':>5} {'h':>10} "
                   f"{'L2(psi)':>15} {'Linf(psi)':>15} "
                   f"{'p_L2':>10} {'p_Linf':>10}\n")
        fout.write("-" * 100 + "\n")

        for (ref_name, nnx, h, L2_psi, Linf_psi), p2, pinf in zip(
                results, rates_psi_L2, rates_psi_Linf):

            p2s   = "   -" if p2   is None else f"{p2:10.3f}"
            pinfs = "   -" if pinf is None else f"{pinf:10.3f}"

            fout.write(
                f"{ref_name:<6} {nnx:5d} {h:10.6f} "
                f"{L2_psi:15.6e} {Linf_psi:15.6e} "
                f"{p2s} {pinfs}\n"
            )

print(f"Wrote convergence table to: {OUT_TXT}")


