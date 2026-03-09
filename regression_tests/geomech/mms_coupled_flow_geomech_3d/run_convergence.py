"""Automated convergence study for MMS coupled flow-geomechanics.

Runs PFLOTRAN at multiple grid refinement levels and computes L2/Linf errors
for both flow (pressure) and geomechanics (displacement).

Usage:
  python run_convergence.py [--levels 4,8,16] [--keep] [--pflotran /path/to/pflotran]
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

# Elastic properties
E = 1.0e9
NU = 0.25
RHO_S = 2500.0
RHO_F = 997.16
PHI = 0.2
BETA = 1.0  # Biot coefficient

LAM = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
MU = E / (2.0 * (1.0 + NU))

# Flow properties
K_PERM = 1.0e-12
MU_L = 1.0e-3

# MMS amplitudes
AU = 1.0
AP = 1000.0
P0 = 200000.0

# Default PFLOTRAN path
PFLOTRAN_DEFAULT = str(
    Path(__file__).resolve().parent.parent.parent.parent / "src" / "pflotran" / "pflotran"
)


# Import the data generation functions from the companion script
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_mms_data import generate_body_force, generate_flow_source


def analytical_displacement(x, y, z):
    """Analytical displacement field at vertex coordinates.

    u = Au * [f, f, f] where f = sin(pi*x)*sin(pi*y)*sin(pi*z).
    """
    f = AU * np.sin(PI * x) * np.sin(PI * y) * np.sin(PI * z)
    return f, f, f


def analytical_pressure(xc, yc, zc):
    """Analytical pressure at cell-center coordinates."""
    return P0 + AP * np.sin(PI * xc) * np.sin(PI * yc) * np.sin(PI * zc)


def update_input_nxyz(template_text, nx, ny, nz):
    """Replace NXYZ line in input deck template."""
    lines = template_text.splitlines()
    out = []
    replaced = False
    for line in lines:
        if line.strip().upper().startswith("NXYZ"):
            out.append(f"    NXYZ {nx} {ny} {nz}")
            replaced = True
        else:
            out.append(line)
    if not replaced:
        raise RuntimeError("NXYZ line not found in template input.")
    return "\n".join(out) + "\n"


def run_case(nx, ny, nz, template_text, base_dir, pflotran, keep_outputs=False):
    """Run one PFLOTRAN case in a temporary subdirectory.

    Each grid level runs in its own temp directory to avoid overwriting
    the default data files (mms_body_force.h5, etc.) in the test directory.
    """
    case_tag = f"mms_coupled_flow_geomech_3d_n{nx}"

    # Create a temporary directory for this run
    tmpdir = Path(tempfile.mkdtemp(prefix=f"mms_n{nx}_", dir=base_dir))
    input_path = tmpdir / f"{case_tag}.in"
    stdout_path = tmpdir / f"{case_tag}.stdout"
    stderr_path = tmpdir / f"{case_tag}.stderr"
    flow_h5_path = tmpdir / f"{case_tag}.h5"
    geomech_h5_path = tmpdir / f"{case_tag}-geomech.h5"

    # Write modified input deck
    input_path.write_text(
        update_input_nxyz(template_text, nx, ny, nz), encoding="ascii"
    )

    # Generate MMS data for this grid level (in the temp directory)
    generate_body_force(nx, ny, nz, E, NU, RHO_S, RHO_F, PHI, AU, AP, BETA,
                        str(tmpdir / "mms_body_force.h5"))
    generate_flow_source(nx, ny, nz, K_PERM, MU_L, RHO_F, AP,
                         str(tmpdir / "mms_flow_source.h5"))

    # Run PFLOTRAN
    print(f"  Running n={nx} ...", end=" ", flush=True)
    with stdout_path.open("w", encoding="ascii") as out, \
         stderr_path.open("w", encoding="ascii") as err:
        try:
            subprocess.check_call(
                [pflotran, "-pflotranin", str(input_path)],
                cwd=str(tmpdir), stdout=out, stderr=err, timeout=600
            )
        except subprocess.CalledProcessError:
            print("FAILED")
            print(f"  See {stderr_path} for details.")
            raise
    print("done")

    # ---- Read flow results ----
    p_l2 = p_linf = float("nan")
    if flow_h5_path.exists():
        with h5py.File(flow_h5_path, "r") as h5:
            time_groups = [k for k in h5.keys() if "Time" in k]
            if time_groups:
                tg = sorted(time_groups)[-1]
                p_num = h5[f"{tg}/Liquid_Pressure [Pa]"][...]
                xc = h5[f"{tg}/X_Coordinate [m]"][...]
                yc = h5[f"{tg}/Y_Coordinate [m]"][...]
                zc = h5[f"{tg}/Z_Coordinate [m]"][...]

                p_ana = analytical_pressure(xc, yc, zc)
                p_err = np.abs(p_num - p_ana)
                p_l2 = math.sqrt(np.mean(p_err**2))
                p_linf = float(np.max(p_err))

    # ---- Read geomech results ----
    u_l2 = u_linf = float("nan")
    if geomech_h5_path.exists():
        with h5py.File(geomech_h5_path, "r") as h5:
            time_groups = [k for k in h5.keys() if "Time" in k]
            if time_groups:
                tg = sorted(time_groups)[-1]

                xv = h5["Domain/X"][...]
                yv = h5["Domain/Y"][...]
                zv = h5["Domain/Z"][...]

                ux_num = h5[f"{tg}/displacement_x [m]"][...]
                uy_num = h5[f"{tg}/displacement_y [m]"][...]
                uz_num = h5[f"{tg}/displacement_z [m]"][...]

                ux_ana, uy_ana, uz_ana = analytical_displacement(xv, yv, zv)

                err_x = ux_num - ux_ana
                err_y = uy_num - uy_ana
                err_z = uz_num - uz_ana
                err_vec = np.sqrt(err_x**2 + err_y**2 + err_z**2)

                u_l2 = math.sqrt(np.mean(err_vec**2))
                u_linf = float(np.max(err_vec))

    # Clean up the temporary directory unless keeping outputs
    if not keep_outputs:
        shutil.rmtree(tmpdir)
    else:
        print(f"  Outputs kept in: {tmpdir}")

    return p_l2, p_linf, u_l2, u_linf


def main():
    parser = argparse.ArgumentParser(
        description="Run MMS coupled flow-geomechanics convergence study."
    )
    parser.add_argument("--levels", default="4,8,16",
                        help="Comma-separated NXYZ levels (cubic).")
    parser.add_argument("--keep", action="store_true",
                        help="Keep per-level output files.")
    parser.add_argument("--pflotran", default=PFLOTRAN_DEFAULT,
                        help="Path to pflotran executable.")
    args = parser.parse_args()

    workdir = Path(__file__).resolve().parent
    template_path = workdir / "mms_coupled_flow_geomech_3d.in"
    template_text = template_path.read_text(encoding="ascii")

    levels = [int(v.strip()) for v in args.levels.split(",") if v.strip()]

    print("MMS Coupled Flow-Geomechanics Convergence Study")
    print("=" * 55)

    results = []
    for n in levels:
        p_l2, p_linf, u_l2, u_linf = run_case(
            n, n, n, template_text, workdir, args.pflotran,
            keep_outputs=args.keep
        )
        h = 1.0 / n
        results.append((n, h, p_l2, p_linf, u_l2, u_linf))

    # ---- Write CSV ----
    csv_path = workdir / "mms_convergence.csv"
    with csv_path.open("w", encoding="ascii") as f:
        f.write("n,h,p_l2,p_linf,u_l2,u_linf,p_order,u_order\n")
        prev = None
        for n, h, p_l2, p_linf, u_l2, u_linf in results:
            p_order = ""
            u_order = ""
            if prev is not None:
                _, h0, p_l2_0, _, u_l2_0, _ = prev
                if p_l2 > 0 and p_l2_0 > 0:
                    p_order = f"{math.log(p_l2 / p_l2_0) / math.log(h / h0):.3f}"
                if u_l2 > 0 and u_l2_0 > 0:
                    u_order = f"{math.log(u_l2 / u_l2_0) / math.log(h / h0):.3f}"
            f.write(f"{n},{h:.8e},{p_l2:.8e},{p_linf:.8e},"
                    f"{u_l2:.8e},{u_linf:.8e},{p_order},{u_order}\n")
            prev = (n, h, p_l2, p_linf, u_l2, u_linf)

    # ---- Print results ----
    print()
    print(f"{'n':>4}  {'h':>10}  {'p_L2':>12}  {'p_order':>8}  {'u_L2':>12}  {'u_order':>8}")
    print("-" * 65)
    prev = None
    for n, h, p_l2, p_linf, u_l2, u_linf in results:
        p_ord_str = "---"
        u_ord_str = "---"
        if prev is not None:
            _, h0, p_l2_0, _, u_l2_0, _ = prev
            if p_l2 > 0 and p_l2_0 > 0:
                p_ord_str = f"{math.log(p_l2 / p_l2_0) / math.log(h / h0):.3f}"
            if u_l2 > 0 and u_l2_0 > 0:
                u_ord_str = f"{math.log(u_l2 / u_l2_0) / math.log(h / h0):.3f}"
        print(f"{n:>4}  {h:>10.6f}  {p_l2:>12.4e}  {p_ord_str:>8}  {u_l2:>12.4e}  {u_ord_str:>8}")
        prev = (n, h, p_l2, p_linf, u_l2, u_linf)

    # ---- Plot convergence ----
    h_vals = [r[1] for r in results]
    p_l2_vals = [r[2] for r in results]
    u_l2_vals = [r[4] for r in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Pressure convergence
    ax1.loglog(h_vals, p_l2_vals, "o-", label="Pressure L2 error")
    # Reference slope
    if len(h_vals) >= 2:
        h_ref = np.array(h_vals)
        ax1.loglog(h_ref, p_l2_vals[0] * (h_ref / h_ref[0])**2, "k--",
                   alpha=0.5, label="O(h^2)")
    ax1.invert_xaxis()
    ax1.set_xlabel("h")
    ax1.set_ylabel("L2 error [Pa]")
    ax1.set_title("Flow Pressure Convergence")
    ax1.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax1.legend()

    # Displacement convergence
    ax2.loglog(h_vals, u_l2_vals, "s-", color="C1",
               label="Displacement L2 error")
    if len(h_vals) >= 2:
        ax2.loglog(h_ref, u_l2_vals[0] * (h_ref / h_ref[0])**2, "k--",
                   alpha=0.5, label="O(h^2)")
    ax2.invert_xaxis()
    ax2.set_xlabel("h")
    ax2.set_ylabel("L2 error [m]")
    ax2.set_title("Geomech Displacement Convergence")
    ax2.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax2.legend()

    fig.suptitle("MMS Coupled Flow-Geomechanics Convergence", fontsize=13)
    fig.tight_layout()
    fig.savefig(workdir / "mms_convergence.png", dpi=150)
    print(f"\nConvergence plot: mms_convergence.png")
    print(f"Convergence data: mms_convergence.csv")


if __name__ == "__main__":
    main()
