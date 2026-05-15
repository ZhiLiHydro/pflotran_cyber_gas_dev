#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

INITIAL_VTK_FILE = Path(__file__).with_name("mms_transient_2d-000.vtk")
FINAL_VTK_FILE = Path(__file__).with_name("mms_transient_2d-001.vtk")
GOLD_FILE = Path(__file__).with_name("mms_transient_2d.regression.gold")


def read_vtk_scalars(path: Path) -> dict[str, np.ndarray]:
    scalars: dict[str, np.ndarray] = {}
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("SCALARS"):
            name = line.split()[1]
            i += 2
            data: list[float] = []
            while i < len(lines):
                cur = lines[i].strip()
                if cur.startswith(("SCALARS", "VECTORS", "FIELD", "POINT_DATA", "CELL_DATA")):
                    i -= 1
                    break
                if cur:
                    data.extend(float(x) for x in cur.split())
                i += 1
            scalars[name] = np.asarray(data, dtype=float)
        i += 1
    return scalars


def read_gold(path: Path) -> dict[str, float]:
    vals: dict[str, float] = {}
    with path.open() as f:
        for line in f:
            if ":" in line:
                k, v = line.split(":", 1)
                vals[k.strip()] = float(v.strip())
    return vals


def rel_l2(numerical: np.ndarray, analytical: np.ndarray) -> float:
    num = np.linalg.norm(numerical - analytical)
    den = np.linalg.norm(analytical)
    return float(num / den) if den != 0.0 else float("inf")


def main() -> int:
    if not FINAL_VTK_FILE.exists():
        print(f"VTK not found: {FINAL_VTK_FILE}")
        return 2

    final_scalars = read_vtk_scalars(FINAL_VTK_FILE)
    ux = final_scalars.get("Stokes_Velocity_X")
    uy = final_scalars.get("Stokes_Velocity_Y")
    pres = final_scalars.get("Liquid_Pressure")
    if ux is None or uy is None or pres is None:
        print("Missing required fields in VTK (Stokes_Velocity_X/Y, Liquid_Pressure)")
        return 2

    ncell = ux.size
    nxy = int(round(math.sqrt(ncell)))
    if nxy * nxy != ncell:
        print(f"Unexpected cell count {ncell}")
        return 2

    xs = (np.arange(nxy) + 0.5) / nxy
    ys = (np.arange(nxy) + 0.5) / nxy
    X, Y = np.meshgrid(xs, ys, indexing="xy")

    pi = math.pi
    ux_exact = pi * np.sin(pi * X) * np.cos(pi * Y)
    uy_exact = -pi * np.cos(pi * X) * np.sin(pi * Y)
    p_exact = np.sin(2.0 * pi * X) * np.sin(2.0 * pi * Y)

    ux = ux.reshape(nxy, nxy)
    uy = uy.reshape(nxy, nxy)
    pres = pres.reshape(nxy, nxy)

    rel_u = rel_l2(ux, ux_exact)
    rel_v = rel_l2(uy, uy_exact)
    rel_p = rel_l2(pres, p_exact)

    initial_rel_u = None
    initial_rel_v = None
    if INITIAL_VTK_FILE.exists():
        initial_scalars = read_vtk_scalars(INITIAL_VTK_FILE)
        ux0 = initial_scalars.get("Stokes_Velocity_X")
        uy0 = initial_scalars.get("Stokes_Velocity_Y")
        if ux0 is not None and uy0 is not None and ux0.size == ncell and uy0.size == ncell:
            ux0 = ux0.reshape(nxy, nxy)
            uy0 = uy0.reshape(nxy, nxy)
            initial_rel_u = rel_l2(ux0, ux_exact)
            initial_rel_v = rel_l2(uy0, uy_exact)

    gold = read_gold(GOLD_FILE)
    tol_u = gold.get("rel_tol_u", 1e-6)
    tol_v = gold.get("rel_tol_v", tol_u)
    tol_p = gold.get("rel_tol_p")
    min_improvement = gold.get("min_relaxation_improvement", 0.25)

    ok = True
    if rel_u > tol_u:
        print(f"rel_u {rel_u} > tol {tol_u}")
        ok = False
    if rel_v > tol_v:
        print(f"rel_v {rel_v} > tol {tol_v}")
        ok = False
    if tol_p is not None and rel_p > tol_p:
        print(f"rel_p {rel_p} > tol {tol_p}")
        ok = False
    if initial_rel_u is not None:
        if not (rel_u < min_improvement * initial_rel_u):
            print(
                f"final rel_u {rel_u} did not improve enough from initial {initial_rel_u}; "
                f"expected factor < {min_improvement}"
            )
            ok = False
    if initial_rel_v is not None:
        if not (rel_v < min_improvement * initial_rel_v):
            print(
                f"final rel_v {rel_v} did not improve enough from initial {initial_rel_v}; "
                f"expected factor < {min_improvement}"
            )
            ok = False

    if initial_rel_u is not None and initial_rel_v is not None:
        print(
            f"initial_rel_u: {initial_rel_u:.6e} initial_rel_v: {initial_rel_v:.6e} "
            f"final_rel_u: {rel_u:.6e} final_rel_v: {rel_v:.6e} rel_p: {rel_p:.6e}"
        )
    else:
        print(f"final_rel_u: {rel_u:.6e} final_rel_v: {rel_v:.6e} rel_p: {rel_p:.6e}")
    return 0 if ok else 3


if __name__ == "__main__":
    raise SystemExit(main())
