"""Compare GEOMECHANICS_FLOW_INTERPOLATION_ORDER 0 vs 1 for MMS Biot 3D.

Runs PFLOTRAN twice (coupling_order=0 and coupling_order=1) on a fine mesh,
then plots midplane contours and diagonal line profiles comparing both
solutions against the analytical MMS field.

Usage:
  python compare_coupling_order_biot.py [--n 32] [--pflotran /path/to/pflotran]
"""

import argparse
import math
import re
import shutil
import subprocess
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

E = 1.0e9
NU = 0.25
RHO = 1000.0
ALPHA = 1.0

RHO_W = 1000.0
G_Z = -9.80665

LAM = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
MU = E / (2.0 * (1.0 + NU))

PFLOTRAN_DEFAULT = str(
    Path(__file__).resolve().parent.parent.parent.parent
    / "src" / "pflotran" / "pflotran"
)


# ---- helpers ---------------------------------------------------------------

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
        sinz, cosz = math.sin(PI * z), math.cos(PI * z)
        for iy in range(nyv):
            y = yv[iy]
            siny, cosy = math.sin(PI * y), math.cos(PI * y)
            for ix in range(nxv):
                x = xv[ix]
                sinx, cosx = math.sin(PI * x), math.cos(PI * x)

                f = sinx * siny * sinz
                lap_f = -3.0 * PI**2 * f
                ddiv_dx = PI**2 * (-f + cosx * cosy * sinz + cosx * siny * cosz)
                ddiv_dy = PI**2 * (-f + cosx * cosy * sinz + sinx * cosy * cosz)
                ddiv_dz = PI**2 * (-f + cosx * siny * cosz + sinx * cosy * cosz)

                body_x[idx] = (-(LAM + MU) * ddiv_dx - MU * lap_f) / RHO
                body_y[idx] = (-(LAM + MU) * ddiv_dy - MU * lap_f) / RHO
                body_z[idx] = (-(LAM + MU) * ddiv_dz - MU * lap_f
                               + ALPHA * dp_dz) / RHO
                idx += 1

    cell_ids = np.arange(1, nnode + 1, dtype=np.int32)
    with h5py.File(outpath, "w") as h5:
        h5.create_dataset("Cell Ids", data=cell_ids)
        h5.create_dataset("Body_force_x", data=body_x)
        h5.create_dataset("Body_force_y", data=body_y)
        h5.create_dataset("Body_force_z", data=body_z)


def patch_input(template_text, nx, coupling_order):
    """Return input-deck text with NXYZ and COUPLING_ORDER replaced."""
    lines = template_text.splitlines()
    out = []
    for line in lines:
        stripped = line.strip().upper()
        if stripped.startswith("NXYZ"):
            out.append(f"    NXYZ {nx} {nx} {nx}")
        elif stripped.startswith("GEOMECHANICS_FLOW_INTERPOLATION_ORDER"):
            out.append(f"  GEOMECHANICS_FLOW_INTERPOLATION_ORDER {coupling_order}")
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def run_case(nx, coupling_order, template_text, pflotran, workdir):
    """Run PFLOTRAN and return (tmpdir, geomech_h5_path)."""
    tag = f"mms_biot_ord{coupling_order}_n{nx}"
    tmpdir = Path(tempfile.mkdtemp(prefix=f"{tag}_", dir=workdir))
    inp = tmpdir / f"{tag}.in"
    inp.write_text(patch_input(template_text, nx, coupling_order),
                   encoding="ascii")
    generate_body_force(nx, nx, nx, str(tmpdir / "mms_body_force.h5"))

    print(f"  Running order={coupling_order}, n={nx} ...", end=" ", flush=True)
    with (tmpdir / f"{tag}.stdout").open("w") as out, \
         (tmpdir / f"{tag}.stderr").open("w") as err:
        subprocess.check_call(
            [pflotran, "-pflotranin", str(inp)],
            cwd=str(tmpdir), stdout=out, stderr=err, timeout=1200,
        )
    print("done")
    return tmpdir, tmpdir / f"{tag}-geomech.h5"


def read_geomech(h5_path):
    """Return (x, y, z, ux, uy, uz) arrays from a geomech HDF5 file."""
    with h5py.File(h5_path, "r") as h5:
        tg = sorted(k for k in h5.keys() if "Time" in k)[-1]
        x = h5["Domain/X"][...]
        y = h5["Domain/Y"][...]
        z = h5["Domain/Z"][...]
        ux = h5[f"{tg}/displacement_x [m]"][...]
        uy = h5[f"{tg}/displacement_y [m]"][...]
        uz = h5[f"{tg}/displacement_z [m]"][...]
    return x, y, z, ux, uy, uz


# ---- plotting --------------------------------------------------------------

def make_plots(xv, yv, zv, data_ord0, data_ord1, nx, workdir):
    """Create comparison figures.

    data_ord0/1 = (ux, uy, uz) numerical displacement arrays.
    """
    ux0, uy0, uz0 = data_ord0
    ux1, uy1, uz1 = data_ord1

    # Analytical
    ux_ana = np.sin(PI * xv) * np.sin(PI * yv) * np.sin(PI * zv)

    # ---- midplane mask ----
    z_mid = 0.5
    tol = 0.5 / nx
    mask = np.abs(zv - z_mid) < tol
    xs, ys = xv[mask], yv[mask]
    tri = Triangulation(xs, ys)

    ana_s = ux_ana[mask]
    o0_s = ux0[mask]
    o1_s = ux1[mask]
    err0_s = o0_s - ana_s
    err1_s = o1_s - ana_s

    # =================================================================
    # Figure 1: Contour comparison (4 columns x 2 rows)
    #   Row 0: u_x contours  — Analytical | Order 0 | Order 1 | Diff (0-1)
    #   Row 1: error contours — Error Ord0 | Error Ord1 | Diff (0-1) zoom | scatter
    # =================================================================
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))

    # Shared color range for the field
    vmin = min(np.min(ana_s), np.min(o0_s), np.min(o1_s))
    vmax = max(np.max(ana_s), np.max(o0_s), np.max(o1_s))

    # Row 0, Col 0: Analytical
    ax = axes[0, 0]
    c = ax.tricontourf(tri, ana_s, levels=32, cmap="RdYlBu_r",
                       vmin=vmin, vmax=vmax)
    fig.colorbar(c, ax=ax, label="u_x [m]")
    ax.set_title("Analytical u_x"); ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_aspect("equal")

    # Row 0, Col 1: Order 0
    ax = axes[0, 1]
    c = ax.tricontourf(tri, o0_s, levels=32, cmap="RdYlBu_r",
                       vmin=vmin, vmax=vmax)
    fig.colorbar(c, ax=ax, label="u_x [m]")
    ax.set_title("Numerical u_x (order 0)"); ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_aspect("equal")

    # Row 0, Col 2: Order 1
    ax = axes[0, 2]
    c = ax.tricontourf(tri, o1_s, levels=32, cmap="RdYlBu_r",
                       vmin=vmin, vmax=vmax)
    fig.colorbar(c, ax=ax, label="u_x [m]")
    ax.set_title("Numerical u_x (order 1)"); ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_aspect("equal")

    # Row 0, Col 3: Difference (order 0 - order 1)
    ax = axes[0, 3]
    diff_s = o0_s - o1_s
    dmax = max(abs(np.min(diff_s)), abs(np.max(diff_s)))
    if dmax == 0:
        dmax = 1e-16
    c = ax.tricontourf(tri, diff_s, levels=32, cmap="RdBu_r",
                       vmin=-dmax, vmax=dmax)
    fig.colorbar(c, ax=ax, label="diff [m]", format="%.2e")
    ax.set_title(f"Order 0 - Order 1\nmax|diff| = {dmax:.3e}")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_aspect("equal")

    # Row 1, Col 0: Error order 0
    ax = axes[1, 0]
    emax0 = max(abs(np.min(err0_s)), abs(np.max(err0_s)))
    if emax0 == 0:
        emax0 = 1e-16
    # Use shared error scale for order 0 and order 1
    emax_shared = max(emax0, max(abs(np.min(err1_s)), abs(np.max(err1_s))))
    c = ax.tricontourf(tri, err0_s, levels=32, cmap="RdBu_r",
                       vmin=-emax_shared, vmax=emax_shared)
    fig.colorbar(c, ax=ax, label="error [m]", format="%.2e")
    ax.set_title(f"Error (order 0)\nLinf = {np.max(np.abs(err0_s)):.3e}")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_aspect("equal")

    # Row 1, Col 1: Error order 1
    ax = axes[1, 1]
    c = ax.tricontourf(tri, err1_s, levels=32, cmap="RdBu_r",
                       vmin=-emax_shared, vmax=emax_shared)
    fig.colorbar(c, ax=ax, label="error [m]", format="%.2e")
    ax.set_title(f"Error (order 1)\nLinf = {np.max(np.abs(err1_s)):.3e}")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_aspect("equal")

    # Row 1, Col 2: Error difference contour (zoomed)
    ax = axes[1, 2]
    err_diff = err0_s - err1_s  # = (o0 - ana) - (o1 - ana) = o0 - o1
    c = ax.tricontourf(tri, err_diff, levels=32, cmap="RdBu_r",
                       vmin=-dmax, vmax=dmax)
    fig.colorbar(c, ax=ax, label="diff [m]", format="%.2e")
    ax.set_title(f"Error(ord0) - Error(ord1)\nmax|diff| = {dmax:.3e}")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_aspect("equal")

    # Row 1, Col 3: Scatter — numerical vs analytical for both orders
    ax = axes[1, 3]
    ax.scatter(ana_s, o0_s, s=6, alpha=0.4, label="Order 0", zorder=2)
    ax.scatter(ana_s, o1_s, s=6, alpha=0.4, label="Order 1", zorder=3)
    lims = [vmin, vmax]
    ax.plot(lims, lims, "k--", linewidth=1, zorder=1)
    ax.set_xlabel("Analytical u_x [m]"); ax.set_ylabel("Numerical u_x [m]")
    ax.set_title("Numerical vs Analytical"); ax.legend(markerscale=3)
    ax.set_aspect("equal")

    fig.suptitle(
        f"MMS Biot 3D — Coupling Order 0 vs 1 — Midplane z = 0.5 — n = {nx}",
        fontsize=14, y=1.01,
    )
    fig.tight_layout()
    p1 = workdir / f"mms_biot_coupling_order_contours_n{nx}.png"
    fig.savefig(p1, dpi=200, bbox_inches="tight")
    print(f"\nContour plot saved: {p1}")

    # =================================================================
    # Figure 2: 1-D diagonal line plots
    # =================================================================
    diag_tol = 1.5 / nx
    mask_d = (np.abs(xv - yv) < diag_tol) & (np.abs(yv - zv) < diag_tol)
    if np.sum(mask_d) < 4:
        return

    xd = xv[mask_d]
    si = np.argsort(xd)
    sd = xd[si]
    ana_d = ux_ana[mask_d][si]
    o0_d = ux0[mask_d][si]
    o1_d = ux1[mask_d][si]

    # Fine reference
    s_ref = np.linspace(0, 1, 500)
    f_ref = np.sin(PI * s_ref) ** 3

    fig2, axes2 = plt.subplots(1, 3, figsize=(17, 5))

    # Panel 1: u_x along diagonal
    ax = axes2[0]
    ax.plot(s_ref, f_ref, "k-", lw=1.5, label="Analytical")
    ax.plot(sd, o0_d, "o", ms=3, alpha=0.6, label="Order 0")
    ax.plot(sd, o1_d, "s", ms=3, alpha=0.6, label="Order 1")
    ax.set_xlabel("x = y = z"); ax.set_ylabel("u_x [m]")
    ax.set_title("u_x along body diagonal"); ax.legend(); ax.grid(True, ls="--", lw=0.4)

    # Panel 2: Error along diagonal
    ax = axes2[1]
    ax.plot(sd, o0_d - ana_d, "o-", ms=3, alpha=0.7, label="Order 0 error")
    ax.plot(sd, o1_d - ana_d, "s-", ms=3, alpha=0.7, label="Order 1 error")
    ax.set_xlabel("x = y = z"); ax.set_ylabel("error [m]")
    ax.set_title("u_x error along diagonal"); ax.legend()
    ax.grid(True, ls="--", lw=0.4)
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(-3, 3))

    # Panel 3: Difference (order 0 - order 1) along diagonal
    ax = axes2[2]
    ax.plot(sd, o0_d - o1_d, "d-", ms=3, color="C2")
    ax.set_xlabel("x = y = z"); ax.set_ylabel("order 0 - order 1 [m]")
    ax.set_title("Difference along diagonal")
    ax.grid(True, ls="--", lw=0.4)
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(-3, 3))

    fig2.suptitle(
        f"MMS Biot 3D — Coupling Order 0 vs 1 — Body Diagonal — n = {nx}",
        fontsize=13, y=1.01,
    )
    fig2.tight_layout()
    p2 = workdir / f"mms_biot_coupling_order_diagonal_n{nx}.png"
    fig2.savefig(p2, dpi=200, bbox_inches="tight")
    print(f"Diagonal plot saved:  {p2}")


# ---- main ------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare coupling order 0 vs 1 for MMS Biot 3D."
    )
    parser.add_argument("--n", type=int, default=32,
                        help="Grid level (default 32)")
    parser.add_argument("--pflotran", default=PFLOTRAN_DEFAULT)
    parser.add_argument("--keep", action="store_true",
                        help="Keep temporary output directories.")
    args = parser.parse_args()

    workdir = Path(__file__).resolve().parent
    template_text = (workdir / "mms_biot_3d.in").read_text(encoding="ascii")

    print(f"MMS Biot 3D — Coupling Order 0 vs 1 (n = {args.n})")
    print("=" * 55)

    # Run both coupling orders
    tmp0, h5_0 = run_case(args.n, 0, template_text, args.pflotran, workdir)
    tmp1, h5_1 = run_case(args.n, 1, template_text, args.pflotran, workdir)

    # Read results
    xv, yv, zv, ux0, uy0, uz0 = read_geomech(h5_0)
    _,  _,  _, ux1, uy1, uz1 = read_geomech(h5_1)

    # Analytical
    ux_ana = np.sin(PI * xv) * np.sin(PI * yv) * np.sin(PI * zv)

    # Print summary
    for label, ux_num in [("Order 0", ux0), ("Order 1", ux1)]:
        err = ux_num - ux_ana
        l2 = math.sqrt(np.mean(err**2))
        linf = float(np.max(np.abs(err)))
        print(f"  {label}:  u_x L2 = {l2:.4e}   Linf = {linf:.4e}")
    diff = ux0 - ux1
    print(f"  Diff (0-1): max|diff| = {np.max(np.abs(diff)):.4e}")

    # Plot
    make_plots(xv, yv, zv, (ux0, uy0, uz0), (ux1, uy1, uz1), args.n, workdir)

    # Cleanup
    if not args.keep:
        shutil.rmtree(tmp0)
        shutil.rmtree(tmp1)
    else:
        print(f"\nOutputs: {tmp0}\n         {tmp1}")


if __name__ == "__main__":
    main()
