import argparse
import math
import os
import subprocess
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

PI = math.pi

E = 1.0e9
NU = 0.25
RHO = 1000.0

LAM = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
MU = E / (2.0 * (1.0 + NU))

PFLOTRAN = "/Users/karr502/src/pflotran-geomech-refactor-part2/src/pflotran/pflotran"


def generate_body_force(nx, ny, nz, filename):
    nxv = nx + 1
    nyv = ny + 1
    nzv = nz + 1
    nxyv = nxv * nyv
    nnode = nxv * nyv * nzv

    xv = np.linspace(0.0, 1.0, nxv)
    yv = np.linspace(0.0, 1.0, nyv)
    zv = np.linspace(0.0, 1.0, nzv)

    body_x = np.zeros(nnode, dtype=np.float64)
    body_y = np.zeros(nnode, dtype=np.float64)
    body_z = np.zeros(nnode, dtype=np.float64)

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
                body_z[idx] = (-(LAM + MU) * ddiv_dz - MU * lap_f) / RHO
                idx += 1

    cell_ids = np.arange(1, nnode + 1, dtype=np.int32)
    with h5py.File(filename, "w") as h5:
        h5.create_dataset("Cell Ids", data=cell_ids)
        h5.create_dataset("Body_force_x", data=body_x)
        h5.create_dataset("Body_force_y", data=body_y)
        h5.create_dataset("Body_force_z", data=body_z)


def update_input_nxyz(template_text, nx, ny, nz):
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


def run_case(nx, ny, nz, template_text, workdir, keep_outputs=False):
    case_tag = f"mms_elastic_3d_n{nx}"
    input_path = workdir / f"{case_tag}.in"
    stdout_path = workdir / f"{case_tag}.stdout"
    stderr_path = workdir / f"{case_tag}.stderr"
    h5_path = workdir / f"{case_tag}-geomech.h5"
    xmf_path = workdir / f"{case_tag}-geomech-000.xmf"

    input_path.write_text(update_input_nxyz(template_text, nx, ny, nz), encoding="ascii")

    generate_body_force(nx, ny, nz, workdir / "mms_body_force.h5")

    with stdout_path.open("w", encoding="ascii") as out, stderr_path.open("w", encoding="ascii") as err:
        subprocess.check_call([PFLOTRAN, "-pflotranin", str(input_path)], cwd=workdir, stdout=out, stderr=err)

    if not h5_path.exists():
        raise RuntimeError(f"Missing geomech output: {h5_path}")

    with h5py.File(h5_path, "r") as h5:
        time_groups = [k for k in h5.keys() if "Time" in k]
        if not time_groups:
            raise RuntimeError("No time groups found in geomech output.")
        time_group = sorted(time_groups)[-1]

        x = h5["Domain/X"][...]
        y = h5["Domain/Y"][...]
        z = h5["Domain/Z"][...]

        ux = h5[f"{time_group}/displacement_x [m]"][...]
        uy = h5[f"{time_group}/displacement_y [m]"][...]
        uz = h5[f"{time_group}/displacement_z [m]"][...]

    ana = np.sin(PI * x) * np.sin(PI * y) * np.sin(PI * z)

    err_x = ux - ana
    err_y = uy - ana
    err_z = uz - ana

    err_vec = np.sqrt(err_x ** 2 + err_y ** 2 + err_z ** 2)
    ana_vec = np.sqrt(ana ** 2 + ana ** 2 + ana ** 2)

    l2 = math.sqrt(np.mean(err_vec ** 2))
    l2_rel = math.sqrt(np.sum(err_vec ** 2) / np.sum(ana_vec ** 2))
    linf = float(np.max(np.abs(err_vec)))

    if not keep_outputs:
        for p in [stdout_path, stderr_path, h5_path, xmf_path]:
            if p.exists():
                p.unlink()
        if input_path.exists():
            input_path.unlink()

    return l2, l2_rel, linf


def main():
    parser = argparse.ArgumentParser(description="Run MMS geomechanics convergence study.")
    parser.add_argument("--levels", default="4,8,16", help="Comma-separated NXYZ levels (cubic).")
    parser.add_argument("--keep", action="store_true", help="Keep per-level output files.")
    args = parser.parse_args()

    workdir = Path(__file__).resolve().parent
    template_path = workdir / "mms_elastic_3d.in"
    template_text = template_path.read_text(encoding="ascii")

    levels = [int(v.strip()) for v in args.levels.split(",") if v.strip()]

    results = []
    for n in levels:
        l2, l2_rel, linf = run_case(n, n, n, template_text, workdir, keep_outputs=args.keep)
        h = 1.0 / n
        results.append((n, h, l2, l2_rel, linf))

    csv_path = workdir / "mms_convergence.csv"
    with csv_path.open("w", encoding="ascii") as f:
        f.write("n,h,l2,l2_rel,linf,order\n")
        prev = None
        for n, h, l2, l2_rel, linf in results:
            order = ""
            if prev is not None:
                _, h0, l20, _, _ = prev
                order_val = math.log(l2 / l20) / math.log(h / h0)
                order = f"{order_val:.3f}"
            f.write(f"{n},{h:.8e},{l2:.8e},{l2_rel:.8e},{linf:.8e},{order}\n")
            prev = (n, h, l2, l2_rel, linf)

    h_vals = [r[1] for r in results]
    l2_vals = [r[2] for r in results]

    plt.figure(figsize=(5, 4))
    plt.loglog(h_vals, l2_vals, "o-", label="L2 error")
    plt.gca().invert_xaxis()
    plt.xlabel("h")
    plt.ylabel("L2 error")
    plt.title("MMS Elastic 3D Convergence")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(workdir / "mms_convergence.png", dpi=150)


if __name__ == "__main__":
    main()
