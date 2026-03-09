"""Compare GEOMECHANICS_FLOW_INTERPOLATION_ORDER 0 vs 1 for MMS coupled flow-geomechanics.

Runs convergence studies for both coupling orders and produces comparison plots.

Usage:
  python compare_coupling_order.py [--levels 4,8,16] [--pflotran /path/to/pflotran]
"""

import argparse
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- MMS parameters (must match input deck and make_mms_data.py) ----
PI = math.pi

E = 1.0e9
NU = 0.25
RHO_S = 2500.0
RHO_F = 997.16
PHI = 0.2
BETA = 1.0

LAM = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
MU_E = E / (2.0 * (1.0 + NU))

K_PERM = 1.0e-12
MU_L = 1.0e-3

AU = 1.0
AP = 1000.0
P0 = 200000.0

PFLOTRAN_DEFAULT = str(
    Path(__file__).resolve().parent.parent.parent.parent / "src" / "pflotran" / "pflotran"
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_mms_data import generate_body_force, generate_flow_source


def analytical_displacement(x, y, z):
    f = AU * np.sin(PI * x) * np.sin(PI * y) * np.sin(PI * z)
    return f, f, f


def analytical_pressure(xc, yc, zc):
    return P0 + AP * np.sin(PI * xc) * np.sin(PI * yc) * np.sin(PI * zc)


def update_input(template_text, nx, ny, nz, coupling_order):
    """Replace NXYZ and GEOMECHANICS_FLOW_INTERPOLATION_ORDER in input deck."""
    lines = template_text.splitlines()
    out = []
    for line in lines:
        stripped = line.strip().upper()
        if stripped.startswith("NXYZ"):
            out.append(f"    NXYZ {nx} {ny} {nz}")
        elif stripped.startswith("GEOMECHANICS_FLOW_INTERPOLATION_ORDER"):
            out.append(f"  GEOMECHANICS_FLOW_INTERPOLATION_ORDER {coupling_order}")
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def run_case(nx, coupling_order, template_text, base_dir, pflotran):
    """Run one PFLOTRAN case, return error metrics."""
    case_tag = f"mms_coupled_flow_geomech_3d_n{nx}"
    tmpdir = Path(tempfile.mkdtemp(prefix=f"mms_co{coupling_order}_n{nx}_", dir=base_dir))

    input_path = tmpdir / f"{case_tag}.in"
    flow_h5_path = tmpdir / f"{case_tag}.h5"
    geomech_h5_path = tmpdir / f"{case_tag}-geomech.h5"

    input_path.write_text(
        update_input(template_text, nx, nx, nx, coupling_order), encoding="ascii"
    )

    generate_body_force(nx, nx, nx, E, NU, RHO_S, RHO_F, PHI, AU, AP, BETA,
                        str(tmpdir / "mms_body_force.h5"))
    generate_flow_source(nx, nx, nx, K_PERM, MU_L, RHO_F, AP,
                         str(tmpdir / "mms_flow_source.h5"))

    stdout_path = tmpdir / "stdout.txt"
    stderr_path = tmpdir / "stderr.txt"

    with stdout_path.open("w") as out, stderr_path.open("w") as err:
        try:
            subprocess.check_call(
                [pflotran, "-pflotranin", str(input_path)],
                cwd=str(tmpdir), stdout=out, stderr=err, timeout=600
            )
        except subprocess.CalledProcessError:
            print(f"FAILED (order={coupling_order}, n={nx})")
            print(f"  See {stderr_path}")
            raise

    # Read flow results
    p_l2 = p_linf = float("nan")
    if flow_h5_path.exists():
        with h5py.File(flow_h5_path, "r") as h5:
            time_groups = sorted([k for k in h5.keys() if "Time" in k])
            if time_groups:
                tg = time_groups[-1]
                p_num = h5[f"{tg}/Liquid_Pressure [Pa]"][...]
                xc = h5[f"{tg}/X_Coordinate [m]"][...]
                yc = h5[f"{tg}/Y_Coordinate [m]"][...]
                zc = h5[f"{tg}/Z_Coordinate [m]"][...]
                p_ana = analytical_pressure(xc, yc, zc)
                p_err = np.abs(p_num - p_ana)
                p_l2 = math.sqrt(np.mean(p_err**2))
                p_linf = float(np.max(p_err))

    # Read geomech results
    u_l2 = u_linf = float("nan")
    if geomech_h5_path.exists():
        with h5py.File(geomech_h5_path, "r") as h5:
            time_groups = sorted([k for k in h5.keys() if "Time" in k])
            if time_groups:
                tg = time_groups[-1]
                xv = h5["Domain/X"][...]
                yv = h5["Domain/Y"][...]
                zv = h5["Domain/Z"][...]
                ux_num = h5[f"{tg}/displacement_x [m]"][...]
                uy_num = h5[f"{tg}/displacement_y [m]"][...]
                uz_num = h5[f"{tg}/displacement_z [m]"][...]
                ux_ana, uy_ana, uz_ana = analytical_displacement(xv, yv, zv)
                err_vec = np.sqrt((ux_num - ux_ana)**2 +
                                  (uy_num - uy_ana)**2 +
                                  (uz_num - uz_ana)**2)
                u_l2 = math.sqrt(np.mean(err_vec**2))
                u_linf = float(np.max(err_vec))

    shutil.rmtree(tmpdir)
    return p_l2, p_linf, u_l2, u_linf


def compute_orders(results):
    """Compute convergence orders from list of (n, h, p_l2, p_linf, u_l2, u_linf)."""
    p_orders = []
    u_orders = []
    for i in range(1, len(results)):
        _, h0, p0, _, u0, _ = results[i - 1]
        _, h1, p1, _, u1, _ = results[i]
        if p0 > 0 and p1 > 0:
            p_orders.append(math.log(p1 / p0) / math.log(h1 / h0))
        else:
            p_orders.append(float("nan"))
        if u0 > 0 and u1 > 0:
            u_orders.append(math.log(u1 / u0) / math.log(h1 / h0))
        else:
            u_orders.append(float("nan"))
    return p_orders, u_orders


def main():
    parser = argparse.ArgumentParser(
        description="Compare coupling_order 0 vs 1 for MMS coupled flow-geomechanics."
    )
    parser.add_argument("--levels", default="4,8,16",
                        help="Comma-separated NXYZ levels (cubic).")
    parser.add_argument("--pflotran", default=PFLOTRAN_DEFAULT,
                        help="Path to pflotran executable.")
    args = parser.parse_args()

    workdir = Path(__file__).resolve().parent
    template_path = workdir / "mms_coupled_flow_geomech_3d.in"
    template_text = template_path.read_text(encoding="ascii")
    levels = [int(v.strip()) for v in args.levels.split(",") if v.strip()]

    all_results = {}

    for co in [0, 1]:
        print(f"\n{'='*60}")
        print(f"  GEOMECHANICS_FLOW_INTERPOLATION_ORDER = {co}")
        print(f"{'='*60}")
        results = []
        for n in levels:
            print(f"  n={n} ...", end=" ", flush=True)
            p_l2, p_linf, u_l2, u_linf = run_case(n, co, template_text, workdir, args.pflotran)
            h = 1.0 / n
            results.append((n, h, p_l2, p_linf, u_l2, u_linf))
            print(f"p_L2={p_l2:.4e}  u_L2={u_l2:.4e}")
        all_results[co] = results

    # ---- Print comparison table ----
    print(f"\n{'='*80}")
    print("COMPARISON TABLE")
    print(f"{'='*80}")

    for co in [0, 1]:
        results = all_results[co]
        p_orders, u_orders = compute_orders(results)
        print(f"\n  Coupling Order = {co}")
        print(f"  {'n':>4}  {'h':>10}  {'p_L2':>12}  {'p_order':>8}  {'u_L2':>12}  {'u_order':>8}")
        print(f"  {'-'*62}")
        for i, (n, h, p_l2, p_linf, u_l2, u_linf) in enumerate(results):
            p_o = f"{p_orders[i-1]:.3f}" if i > 0 else "---"
            u_o = f"{u_orders[i-1]:.3f}" if i > 0 else "---"
            print(f"  {n:>4}  {h:>10.6f}  {p_l2:>12.4e}  {p_o:>8}  {u_l2:>12.4e}  {u_o:>8}")

    # ---- Generate plots ----
    h_vals = [r[1] for r in all_results[0]]
    h_arr = np.array(h_vals)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    colors = {0: "C0", 1: "C1"}
    markers = {0: "o", 1: "s"}
    labels = {0: "Coupling Order 0 (flow first)",
              1: "Coupling Order 1 (geomech first)"}

    # ---- Plot 1: Pressure L2 error ----
    ax = axes[0, 0]
    for co in [0, 1]:
        p_l2_vals = [r[2] for r in all_results[co]]
        ax.loglog(h_vals, p_l2_vals, f"{markers[co]}-", color=colors[co],
                  label=labels[co], markersize=8, linewidth=2)
    # Reference O(h^2)
    ref_p = all_results[0][0][2]
    ax.loglog(h_arr, ref_p * (h_arr / h_arr[0])**2, "k--", alpha=0.4,
              linewidth=1.5, label="O(h$^2$)")
    ax.invert_xaxis()
    ax.set_xlabel("h (mesh spacing)", fontsize=12)
    ax.set_ylabel("L$_2$ error [Pa]", fontsize=12)
    ax.set_title("Pressure L$_2$ Error", fontsize=13)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(fontsize=10)

    # ---- Plot 2: Pressure Linf error ----
    ax = axes[0, 1]
    for co in [0, 1]:
        p_linf_vals = [r[3] for r in all_results[co]]
        ax.loglog(h_vals, p_linf_vals, f"{markers[co]}-", color=colors[co],
                  label=labels[co], markersize=8, linewidth=2)
    ref_pinf = all_results[0][0][3]
    ax.loglog(h_arr, ref_pinf * (h_arr / h_arr[0])**2, "k--", alpha=0.4,
              linewidth=1.5, label="O(h$^2$)")
    ax.invert_xaxis()
    ax.set_xlabel("h (mesh spacing)", fontsize=12)
    ax.set_ylabel("L$_\\infty$ error [Pa]", fontsize=12)
    ax.set_title("Pressure L$_\\infty$ Error", fontsize=13)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(fontsize=10)

    # ---- Plot 3: Displacement L2 error ----
    ax = axes[1, 0]
    for co in [0, 1]:
        u_l2_vals = [r[4] for r in all_results[co]]
        ax.loglog(h_vals, u_l2_vals, f"{markers[co]}-", color=colors[co],
                  label=labels[co], markersize=8, linewidth=2)
    ref_u = all_results[0][0][4]
    ax.loglog(h_arr, ref_u * (h_arr / h_arr[0])**2, "k--", alpha=0.4,
              linewidth=1.5, label="O(h$^2$)")
    ax.invert_xaxis()
    ax.set_xlabel("h (mesh spacing)", fontsize=12)
    ax.set_ylabel("L$_2$ error [m]", fontsize=12)
    ax.set_title("Displacement L$_2$ Error", fontsize=13)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(fontsize=10)

    # ---- Plot 4: Displacement Linf error ----
    ax = axes[1, 1]
    for co in [0, 1]:
        u_linf_vals = [r[5] for r in all_results[co]]
        ax.loglog(h_vals, u_linf_vals, f"{markers[co]}-", color=colors[co],
                  label=labels[co], markersize=8, linewidth=2)
    ref_uinf = all_results[0][0][5]
    ax.loglog(h_arr, ref_uinf * (h_arr / h_arr[0])**2, "k--", alpha=0.4,
              linewidth=1.5, label="O(h$^2$)")
    ax.invert_xaxis()
    ax.set_xlabel("h (mesh spacing)", fontsize=12)
    ax.set_ylabel("L$_\\infty$ error [m]", fontsize=12)
    ax.set_title("Displacement L$_\\infty$ Error", fontsize=13)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(fontsize=10)

    fig.suptitle("MMS Coupled Flow-Geomechanics: Coupling Order 0 vs 1",
                 fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    plot_path = workdir / "coupling_order_comparison.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"\nComparison plot saved: {plot_path}")

    # ---- Bar chart of convergence orders ----
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    for co in [0, 1]:
        p_orders, u_orders = compute_orders(all_results[co])
        bar_labels = [f"n={all_results[co][i][0]}->{all_results[co][i+1][0]}"
                      for i in range(len(p_orders))]
        x_pos = np.arange(len(bar_labels))
        width = 0.35
        offset = -width/2 if co == 0 else width/2

        ax1.bar(x_pos + offset, p_orders, width, label=labels[co], color=colors[co], alpha=0.8)
        ax2.bar(x_pos + offset, u_orders, width, label=labels[co], color=colors[co], alpha=0.8)

    ax1.axhline(y=2.0, color="k", linestyle="--", alpha=0.5, linewidth=1.5, label="Ideal (2nd order)")
    ax2.axhline(y=2.0, color="k", linestyle="--", alpha=0.5, linewidth=1.5, label="Ideal (2nd order)")

    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(bar_labels)
    ax1.set_ylabel("Convergence Order", fontsize=12)
    ax1.set_title("Pressure Convergence Order", fontsize=13)
    ax1.legend(fontsize=9)
    ax1.set_ylim(0, 3)
    ax1.grid(axis="y", linestyle="--", alpha=0.5)

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(bar_labels)
    ax2.set_ylabel("Convergence Order", fontsize=12)
    ax2.set_title("Displacement Convergence Order", fontsize=13)
    ax2.legend(fontsize=9)
    ax2.set_ylim(0, 3)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)

    fig2.suptitle("Convergence Orders: Coupling Order 0 vs 1",
                  fontsize=14, fontweight="bold")
    fig2.tight_layout(rect=[0, 0, 1, 0.93])

    orders_path = workdir / "coupling_order_rates.png"
    fig2.savefig(orders_path, dpi=150)
    plt.close(fig2)
    print(f"Convergence orders plot saved: {orders_path}")


if __name__ == "__main__":
    main()
