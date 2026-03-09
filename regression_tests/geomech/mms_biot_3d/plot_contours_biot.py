"""Plot contour comparisons (analytical vs numerical vs error) for MMS Biot 3D.

Runs PFLOTRAN on a fine mesh, then plots midplane (z=0.5) contours for:
  - Displacement (u_x component)
  - Pressure

Each field gets a row of 3 panels: Analytical | Numerical | Error

Usage:
  python plot_contours_biot.py [--n 32] [--pflotran /path/to/pflotran]
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
from matplotlib.tri import Triangulation

# ---- MMS parameters (must match input deck / make_mms_body_force.py) ----
PI = math.pi

# Elastic properties
E = 1.0e9
NU = 0.25
RHO = 1000.0
ALPHA = 1.0

# Hydrostatic pressure
RHO_W = 1000.0
G_Z = -9.80665
P0 = 101325.0
Z0 = 0.0

LAM = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
MU = E / (2.0 * (1.0 + NU))

PFLOTRAN_DEFAULT = str(
    Path(__file__).resolve().parent.parent.parent.parent / "src" / "pflotran" / "pflotran"
)


def generate_body_force(nx, ny, nz, outpath):
    """Generate mms_body_force.h5 for the given grid level."""
    nxv, nyv, nzv = nx + 1, ny + 1, nz + 1
    nnode = nxv * nyv * nzv

    xv = np.linspace(0.0, 1.0, nxv)
    yv = np.linspace(0.0, 1.0, nyv)
    zv = np.linspace(0.0, 1.0, nzv)

    body_x = np.zeros(nnode, dtype=np.float64)
    body_y = np.zeros(nnode, dtype=np.float64)
    body_z = np.zeros(nnode, dtype=np.float64)

    dp_dz = RHO_W * G_Z

    idx = 0
    for iz in range(nzv):
        z = zv[iz]
        sinz = math.sin(PI * z)
        cosz = math.cos(PI * z)
        for iy in range(nyv):
            y = yv[iy]
            siny = math.sin(PI * y)
            cosy = math.cos(PI * y)
            for ix in range(nxv):
                x = xv[ix]
                sinx = math.sin(PI * x)
                cosx = math.cos(PI * x)

                f = sinx * siny * sinz
                lap_f = -3.0 * (PI ** 2) * f

                ddiv_dx = (PI ** 2) * (-f + cosx * cosy * sinz + cosx * siny * cosz)
                ddiv_dy = (PI ** 2) * (-f + cosx * cosy * sinz + sinx * cosy * cosz)
                ddiv_dz = (PI ** 2) * (-f + cosx * siny * cosz + sinx * cosy * cosz)

                body_x[idx] = (-(LAM + MU) * ddiv_dx - MU * lap_f) / RHO
                body_y[idx] = (-(LAM + MU) * ddiv_dy - MU * lap_f) / RHO
                body_z[idx] = (-(LAM + MU) * ddiv_dz - MU * lap_f + ALPHA * dp_dz) / RHO

                idx += 1

    cell_ids = np.arange(1, nnode + 1, dtype=np.int32)
    with h5py.File(outpath, "w") as h5:
        h5.create_dataset("Cell Ids", data=cell_ids)
        h5.create_dataset("Body_force_x", data=body_x)
        h5.create_dataset("Body_force_y", data=body_y)
        h5.create_dataset("Body_force_z", data=body_z)


def update_input_nxyz(template_text, nx, ny, nz):
    """Replace NXYZ line in input deck template."""
    lines = template_text.splitlines()
    out = []
    for line in lines:
        if line.strip().upper().startswith("NXYZ"):
            out.append(f"    NXYZ {nx} {ny} {nz}")
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def run_pflotran(nx, pflotran, workdir):
    """Run PFLOTRAN for the Biot MMS test at a given grid level.

    Returns (geomech_h5_path, flow_h5_path) in the temp directory.
    """
    template_path = workdir / "mms_biot_3d.in"
    template_text = template_path.read_text(encoding="ascii")

    tmpdir = Path(tempfile.mkdtemp(prefix=f"mms_biot_n{nx}_", dir=workdir))
    case_tag = f"mms_biot_3d_n{nx}"
    input_path = tmpdir / f"{case_tag}.in"
    stdout_path = tmpdir / f"{case_tag}.stdout"
    stderr_path = tmpdir / f"{case_tag}.stderr"
    flow_h5 = tmpdir / f"{case_tag}.h5"
    geomech_h5 = tmpdir / f"{case_tag}-geomech.h5"

    input_path.write_text(
        update_input_nxyz(template_text, nx, nx, nx), encoding="ascii"
    )

    generate_body_force(nx, nx, nx, str(tmpdir / "mms_body_force.h5"))

    print(f"  Running PFLOTRAN with n={nx} ...", end=" ", flush=True)
    with stdout_path.open("w") as out, stderr_path.open("w") as err:
        subprocess.check_call(
            [pflotran, "-pflotranin", str(input_path)],
            cwd=str(tmpdir), stdout=out, stderr=err, timeout=1200
        )
    print("done")

    return tmpdir, geomech_h5, flow_h5


def make_contour_plots(tmpdir, geomech_h5, flow_h5, nx, workdir):
    """Generate contour comparison plots from PFLOTRAN output."""

    nxv = nx + 1

    # ---- Read geomech output ----
    with h5py.File(geomech_h5, "r") as h5:
        time_groups = sorted([k for k in h5.keys() if "Time" in k])
        tg = time_groups[-1]
        xv = h5["Domain/X"][...]
        yv = h5["Domain/Y"][...]
        zv = h5["Domain/Z"][...]
        ux_num = h5[f"{tg}/displacement_x [m]"][...]
        uy_num = h5[f"{tg}/displacement_y [m]"][...]
        uz_num = h5[f"{tg}/displacement_z [m]"][...]

    # ---- Analytical displacement ----
    ux_ana = np.sin(PI * xv) * np.sin(PI * yv) * np.sin(PI * zv)
    # u_y = u_z = u_x by construction
    uy_ana = ux_ana.copy()
    uz_ana = ux_ana.copy()

    # Displacement magnitude
    umag_num = np.sqrt(ux_num**2 + uy_num**2 + uz_num**2)
    umag_ana = np.sqrt(ux_ana**2 + uy_ana**2 + uz_ana**2)

    # =====================================================================
    # Midplane slice z ≈ 0.5 (vertex data)
    # =====================================================================
    z_mid = 0.5
    tol = 0.5 / nx  # half an element
    mask_v = np.abs(zv - z_mid) < tol

    xv_s = xv[mask_v]
    yv_s = yv[mask_v]
    tri_v = Triangulation(xv_s, yv_s)

    # Sliced fields
    fields = [
        ("u_x", "displacement u_x", "[m]", "RdYlBu_r",
         ux_ana[mask_v], ux_num[mask_v]),
        ("|u|", "displacement magnitude |u|", "[m]", "magma",
         umag_ana[mask_v], umag_num[mask_v]),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    for row, (tag, label, unit, cmap, ana_s, num_s) in enumerate(fields):
        err_s = num_s - ana_s

        vmin = min(np.min(ana_s), np.min(num_s))
        vmax = max(np.max(ana_s), np.max(num_s))

        # Analytical
        ax = axes[row, 0]
        c = ax.tricontourf(tri_v, ana_s, levels=32, cmap=cmap,
                           vmin=vmin, vmax=vmax)
        fig.colorbar(c, ax=ax, label=f"{tag} {unit}")
        ax.set_title(f"Analytical {label}")
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_aspect("equal")

        # Numerical
        ax = axes[row, 1]
        c = ax.tricontourf(tri_v, num_s, levels=32, cmap=cmap,
                           vmin=vmin, vmax=vmax)
        fig.colorbar(c, ax=ax, label=f"{tag} {unit}")
        ax.set_title(f"Numerical {label}")
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_aspect("equal")

        # Error
        ax = axes[row, 2]
        err_abs_max = max(abs(np.min(err_s)), abs(np.max(err_s)))
        if err_abs_max == 0:
            err_abs_max = 1e-16
        c = ax.tricontourf(tri_v, err_s, levels=32, cmap="RdBu_r",
                           vmin=-err_abs_max, vmax=err_abs_max)
        fig.colorbar(c, ax=ax, label=f"error {unit}", format="%.2e")
        ax.set_title(f"Error (num - ana)\n"
                     f"Linf = {np.max(np.abs(err_s)):.3e}")
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_aspect("equal")

    fig.suptitle(f"MMS Biot 3D — Midplane z = 0.5 — n = {nx}",
                 fontsize=14, y=1.01)
    fig.tight_layout()
    contour_path = workdir / f"mms_biot_3d_contours_n{nx}.png"
    fig.savefig(contour_path, dpi=200, bbox_inches="tight")
    print(f"\nContour plot saved: {contour_path}")

    # =====================================================================
    # 1-D line plots along the body diagonal (0,0,0)→(1,1,1)
    # =====================================================================
    # Extract vertices close to the body diagonal x = y = z
    diag_tol = 1.5 / nx  # ~1.5 element widths
    mask_diag = (np.abs(xv - yv) < diag_tol) & (np.abs(yv - zv) < diag_tol)

    if np.sum(mask_diag) > 3:
        xd = xv[mask_diag]
        # Parameter along diagonal: s = x (since x=y=z on the diagonal)
        sort_idx = np.argsort(xd)
        sd = xd[sort_idx]

        ux_ana_d = ux_ana[mask_diag][sort_idx]
        ux_num_d = ux_num[mask_diag][sort_idx]
        umag_ana_d = umag_ana[mask_diag][sort_idx]
        umag_num_d = umag_num[mask_diag][sort_idx]

        # Fine analytical reference
        s_ref = np.linspace(0, 1, 500)
        f_ref = np.sin(PI * s_ref)**3  # sin(pi*s)^3 on the diagonal
        umag_ref = np.sqrt(3) * np.abs(f_ref)

        fig2, axes2 = plt.subplots(1, 3, figsize=(16, 4.5))

        # u_x along diagonal
        ax = axes2[0]
        ax.plot(s_ref, f_ref, "k-", linewidth=1.5, label="Analytical")
        ax.plot(sd, ux_num_d, "o", markersize=3, alpha=0.7, label="Numerical")
        ax.set_xlabel("x = y = z")
        ax.set_ylabel("u_x [m]")
        ax.set_title("u_x along body diagonal")
        ax.legend()
        ax.grid(True, linestyle="--", linewidth=0.4)

        # |u| along diagonal
        ax = axes2[1]
        ax.plot(s_ref, umag_ref, "k-", linewidth=1.5, label="Analytical")
        ax.plot(sd, umag_num_d, "o", markersize=3, alpha=0.7,
                color="C1", label="Numerical")
        ax.set_xlabel("x = y = z")
        ax.set_ylabel("|u| [m]")
        ax.set_title("|u| along body diagonal")
        ax.legend()
        ax.grid(True, linestyle="--", linewidth=0.4)

        # Error along diagonal
        ax = axes2[2]
        ax.plot(sd, ux_num_d - ux_ana_d, "s-", markersize=3,
                label="u_x error")
        ax.plot(sd, umag_num_d - umag_ana_d, "o-", markersize=3,
                color="C1", label="|u| error")
        ax.set_xlabel("x = y = z")
        ax.set_ylabel("error [m]")
        ax.set_title("Error along body diagonal")
        ax.legend()
        ax.grid(True, linestyle="--", linewidth=0.4)
        ax.ticklabel_format(axis="y", style="scientific", scilimits=(-3, 3))

        fig2.suptitle(f"MMS Biot 3D — Body Diagonal — n = {nx}",
                      fontsize=13, y=1.01)
        fig2.tight_layout()
        line_path = workdir / f"mms_biot_3d_diagonal_n{nx}.png"
        fig2.savefig(line_path, dpi=200, bbox_inches="tight")
        print(f"Diagonal plot saved:  {line_path}")

    # ---- Print global errors ----
    ux_err_all = ux_num - ux_ana
    uy_err_all = uy_num - uy_ana
    uz_err_all = uz_num - uz_ana
    err_vec = np.sqrt(ux_err_all**2 + uy_err_all**2 + uz_err_all**2)

    print(f"\nGlobal displacement errors (n={nx}, h={1.0/nx:.4f}):")
    print(f"  u_x  L2 = {math.sqrt(np.mean(ux_err_all**2)):.4e}   "
          f"Linf = {np.max(np.abs(ux_err_all)):.4e}")
    print(f"  |u|  L2 = {math.sqrt(np.mean(err_vec**2)):.4e}   "
          f"Linf = {np.max(err_vec):.4e}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot MMS Biot 3D contour comparisons (analytical / numerical / error)."
    )
    parser.add_argument("--n", type=int, default=32,
                        help="Grid refinement level (cubic). Default: 32")
    parser.add_argument("--pflotran", default=PFLOTRAN_DEFAULT,
                        help="Path to pflotran executable.")
    parser.add_argument("--keep", action="store_true",
                        help="Keep temporary PFLOTRAN output directory.")
    args = parser.parse_args()

    workdir = Path(__file__).resolve().parent
    print(f"MMS Biot 3D Contour Comparison (n = {args.n})")
    print("=" * 50)

    tmpdir, geomech_h5, flow_h5 = run_pflotran(args.n, args.pflotran, workdir)
    make_contour_plots(tmpdir, geomech_h5, flow_h5, args.n, workdir)

    if not args.keep:
        shutil.rmtree(tmpdir)
    else:
        print(f"Outputs kept in: {tmpdir}")


if __name__ == "__main__":
    main()
