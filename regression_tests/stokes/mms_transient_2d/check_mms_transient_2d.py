#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

VTK_FILE = Path(__file__).with_name("mms_transient_2d-001.vtk")
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
    if not VTK_FILE.exists():
        print(f"VTK not found: {VTK_FILE}")
        return 2

    scalars = read_vtk_scalars(VTK_FILE)
    ux = scalars.get("Stokes_Velocity_X")
    uy = scalars.get("Stokes_Velocity_Y")
    pres = scalars.get("Liquid_Pressure")
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

    gold = read_gold(GOLD_FILE)
    tol_u = gold.get("rel_tol_u", 1e-6)
    tol_v = gold.get("rel_tol_v", tol_u)
    tol_p = gold.get("rel_tol_p")

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

    print(f"rel_u: {rel_u:.6e} rel_v: {rel_v:.6e} rel_p: {rel_p:.6e}")
    return 0 if ok else 3


if __name__ == "__main__":
    raise SystemExit(main())
