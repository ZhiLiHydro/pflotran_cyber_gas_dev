#!/usr/bin/env python3
"""
Post-process check for the Stokes Poiseuille regression test.

Setup: NX=10, NY=10, NZ=1, domain 2.0 x 1.0 x 1.0 m.
       LIQUID_PRESSURE = 1 at west, 0 at east.
       VELOCITY = 0 (no-slip) at south (j=1) and north (j=NY=10).
       VISCOSITY = DENSITY = 1.

Expected: fully-developed pressure-driven flow with a discrete parabolic
          velocity profile.  For no-slip applied at cell centers j=1 and
          j=NY the exact discrete solution is

              u_x[j] = K * (j - 1) * (NY - j)    j = 1 .. NY

          where K depends on the pressure gradient and viscosity.
          In particular u_x[1] = u_x[10] = 0 (no-slip wall cells).

Cell layout (NX=10, NY=10, NZ=1, 1-indexed):
  cell_id = i + (j-1)*NX

Column at i=5 (x = 0.9 m): cells 5, 15, 25, 35, 45, 55, 65, 75, 85, 95
Row at j=5 (y = 0.45 m):   cells 41 .. 50
"""

import math
import re
import sys
from pathlib import Path

# Grid parameters
NX = 10
NY = 10

# Tolerances
PARABOLA_TOL = 1.0e-2   # relative tolerance on parabolic shape ratios
                          # (~0.2% variation seen across j for column at x=0.9;
                          #  allow 1% to be robust against different grid sizes)
UY_MAX = 5.0e-4          # max |u_y| (small non-zero artefact from inlet BCs)
UZ_MAX = 1.0e-12         # u_z must be machine zero

# Cells in the i=5 column (x midpoint, one cell from each j-row)
COLUMN_CELLS = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95]  # j = 1 .. 10

# Cells in the j=5 row (y = 0.45 m)
ROW_CELLS = list(range(41, 51))  # cells 41 .. 50


def parse_regression(path):
    sections = {}
    section = None
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("-- ") and line.endswith(" --"):
                section = line
                sections[section] = {}
                continue
            if section is None:
                continue
            m = re.match(
                r"^(\d+):\s+([+-]?[0-9]*\.?[0-9]+(?:[Ee][+-]?\d+)?)$", line
            )
            if m:
                sections[section][int(m.group(1))] = float(m.group(2))
    return sections


def require(cond, msg):
    if not cond:
        print(f"FAIL: {msg}")
        sys.exit(1)


def main():
    regression_file = (
        Path(__file__).resolve().parent / "stokes_poiseuille.regression"
    )
    sections = parse_regression(str(regression_file))

    p_key  = "-- PRESSURE: Liquid Pressure --"
    ux_key = "-- GENERIC: Stokes Velocity X --"
    uy_key = "-- GENERIC: Stokes Velocity Y --"
    uz_key = "-- GENERIC: Stokes Velocity Z --"

    for key in (p_key, ux_key, uy_key, uz_key):
        require(key in sections, f"Missing section: {key}")

    pressure = sections[p_key]
    vel_x    = sections[ux_key]
    vel_y    = sections[uy_key]
    vel_z    = sections[uz_key]

    # ------------------------------------------------------------------
    # 1. No-slip walls: u_x must be zero at j=1 (cell 5) and j=10 (cell 95)
    # ------------------------------------------------------------------
    require(vel_x[5]  == 0.0, "No-slip violated at south wall (cell 5, j=1).")
    require(vel_x[95] == 0.0, "No-slip violated at north wall (cell 95, j=10).")

    # ------------------------------------------------------------------
    # 2. Parabolic profile in the i=5 column
    #    u_x[j] = K * (j-1) * (NY-j)
    #    Verify using velocity ratios: u_x[j] / u_x[j_ref]  should equal
    #    (j-1)*(NY-j) / (j_ref-1)*(NY-j_ref) to tolerance PARABOLA_TOL.
    #    Use j_ref=5 (near peak) as the reference row.
    # ------------------------------------------------------------------
    j_ref = 5
    cell_ref = COLUMN_CELLS[j_ref - 1]  # cell 45
    ux_ref = vel_x[cell_ref]
    weight_ref = (j_ref - 1) * (NY - j_ref)
    require(ux_ref > 0.0, f"Reference cell {cell_ref} (j={j_ref}) has non-positive u_x.")

    for idx, cell_id in enumerate(COLUMN_CELLS):
        j = idx + 1  # 1-indexed row
        weight = (j - 1) * (NY - j)
        ux = vel_x[cell_id]
        if weight == 0:
            require(
                abs(ux) < 1.0e-14,
                f"Wall cell {cell_id} (j={j}) has non-zero u_x = {ux}.",
            )
        else:
            ratio_sim  = ux / ux_ref
            ratio_anal = weight / weight_ref
            rel_err = abs(ratio_sim - ratio_anal) / ratio_anal
            require(
                rel_err < PARABOLA_TOL,
                f"Parabolic shape violated at j={j} (cell {cell_id}): "
                f"u_x ratio {ratio_sim:.6f} vs analytical {ratio_anal:.6f} "
                f"(rel err {rel_err:.2e}).",
            )

    # ------------------------------------------------------------------
    # 3. Symmetry: u_x[j] == u_x[NY+1-j]  (should be exact by symmetry)
    # ------------------------------------------------------------------
    for j in range(1, NY // 2 + 1):
        cell_lo = COLUMN_CELLS[j - 1]        # j
        cell_hi = COLUMN_CELLS[NY - j]       # NY+1-j
        ux_lo = vel_x[cell_lo]
        ux_hi = vel_x[cell_hi]
        denom = max(abs(ux_lo), abs(ux_hi), 1.0e-30)
        require(
            abs(ux_lo - ux_hi) / denom < 1.0e-10,
            f"Asymmetry in u_x: cell {cell_lo} (j={j}) = {ux_lo:.6e}, "
            f"cell {cell_hi} (j={NY+1-j}) = {ux_hi:.6e}.",
        )

    # ------------------------------------------------------------------
    # 4. Flow direction: all u_x >= 0
    # ------------------------------------------------------------------
    for cell_id, ux in vel_x.items():
        require(ux >= -1.0e-14, f"Negative u_x = {ux:.3e} at cell {cell_id}.")

    # ------------------------------------------------------------------
    # 5. Cross-flow (u_y) should be small
    # ------------------------------------------------------------------
    uy_max_actual = max(abs(v) for v in vel_y.values())
    require(
        uy_max_actual < UY_MAX,
        f"u_y too large: max |u_y| = {uy_max_actual:.3e} (threshold {UY_MAX:.3e}).",
    )

    # ------------------------------------------------------------------
    # 6. u_z must be exactly zero (2-D problem)
    # ------------------------------------------------------------------
    uz_max_actual = max(abs(v) for v in vel_z.values())
    require(
        uz_max_actual < UZ_MAX,
        f"u_z is non-zero: max |u_z| = {uz_max_actual:.3e}.",
    )

    # ------------------------------------------------------------------
    # 7. Pressure monotone decreasing west -> east along j=5 row
    # ------------------------------------------------------------------
    p_row = [pressure[c] for c in ROW_CELLS]
    for idx in range(len(p_row) - 1):
        require(
            p_row[idx] > p_row[idx + 1],
            f"Pressure not monotone at columns {idx+1}->{idx+2} along j=5 row: "
            f"{p_row[idx]:.6e} -> {p_row[idx+1]:.6e}.",
        )

    # ------------------------------------------------------------------
    # 8. u_x nearly uniform along j=5 row (fully-developed flow)
    # ------------------------------------------------------------------
    ux_row = [vel_x[c] for c in ROW_CELLS[1:-1]]  # skip west/east boundary columns
    ux_row_max = max(ux_row)
    ux_row_min = min(ux_row)
    require(
        ux_row_max > 0.0,
        "u_x in j=5 row is zero — no flow detected.",
    )
    variation = (ux_row_max - ux_row_min) / ux_row_max
    require(
        variation < 0.1,
        f"u_x varies too much along j=5 row (variation {variation:.2%}); "
        "flow may not be fully developed.",
    )

    print("stokes_poiseuille checks passed.")


if __name__ == "__main__":
    main()
