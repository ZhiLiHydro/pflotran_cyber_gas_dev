"""
generate_input.py
=================
Generate a PFLOTRAN Stokes input file for flow through a 2D pore network
with randomly packed circular grains.

Grid: 80 x 40 x 1 cells, domain 2 m x 1 m x 0.1 m, dx = dy = 0.025 m.

Grain placement uses Random Sequential Addition (RSA) with a fixed seed.
Grains near the inlet/outlet columns are excluded to guarantee open flow paths.
Connectivity (inlet column -> outlet column) is verified via BFS before writing.
"""

import sys
import numpy as np
from collections import deque

# ---------------------------------------------------------------------------
# Grid parameters
# ---------------------------------------------------------------------------
P_IN = 4.0   # inlet pressure [Pa]; outlet is always 0 Pa

NX, NY, NZ = 80, 40, 1
LX, LY, LZ = 2.0, 1.0, 0.1
DX = LX / NX   # 0.025 m
DY = LY / NY   # 0.025 m
DZ = LZ / NZ   # 0.1 m

# Cell centres
XC = (np.arange(NX) + 0.5) * DX  # shape (NX,)
YC = (np.arange(NY) + 0.5) * DY  # shape (NY,)

# ---------------------------------------------------------------------------
# Grain placement parameters
# ---------------------------------------------------------------------------
SEED      = 7
N_TARGET  = 24       # desired number of grains
R_MIN     = 0.06     # m
R_MAX     = 0.13     # m
CLEARANCE = 0.025    # m  minimum gap between grain surfaces
X_GUARD   = 0.10     # m  grain-free buffer at west and east boundaries

MAX_ATTEMPTS = 100_000   # RSA rejection limit

# ---------------------------------------------------------------------------
# Random Sequential Addition
# ---------------------------------------------------------------------------
rng = np.random.default_rng(seed=SEED)

grains = []   # list of (cx, cy, r)

for _ in range(MAX_ATTEMPTS):
    if len(grains) >= N_TARGET:
        break
    r  = rng.uniform(R_MIN, R_MAX)
    cx = rng.uniform(X_GUARD + r, LX - X_GUARD - r)
    cy = rng.uniform(r, LY - r)

    # Check against domain walls (already ensured by bounds above)
    # Check against all existing grains
    ok = True
    for (ex, ey, er) in grains:
        dist = np.hypot(cx - ex, cy - ey)
        if dist < r + er + CLEARANCE:
            ok = False
            break
    if ok:
        grains.append((cx, cy, r))

print(f"Placed {len(grains)} grains (target {N_TARGET}).")

# ---------------------------------------------------------------------------
# Rasterise grains onto the grid
# ---------------------------------------------------------------------------
# grain_mask[j, i] = True  => cell (i,j) is solid
grain_mask = np.zeros((NY, NX), dtype=bool)

for (cx, cy, r) in grains:
    # Bounding box in cell indices
    i_lo = max(0,    int(np.floor((cx - r) / DX)))
    i_hi = min(NX-1, int(np.floor((cx + r) / DX)))
    j_lo = max(0,    int(np.floor((cy - r) / DY)))
    j_hi = min(NY-1, int(np.floor((cy + r) / DY)))
    for j in range(j_lo, j_hi + 1):
        for i in range(i_lo, i_hi + 1):
            if (XC[i] - cx)**2 + (YC[j] - cy)**2 <= r**2:
                grain_mask[j, i] = True

n_solid = grain_mask.sum()
n_fluid = NX * NY - n_solid
print(f"Solid cells: {n_solid}  Fluid cells: {n_fluid}  "
      f"Porosity: {n_fluid / (NX * NY):.3f}")

# ---------------------------------------------------------------------------
# Connectivity check: BFS from west column to east column
# ---------------------------------------------------------------------------
def bfs_connected(mask):
    """Return True if any fluid cell in column 0 connects to column NX-1."""
    visited = np.zeros((NY, NX), dtype=bool)
    queue   = deque()
    for j in range(NY):
        if not mask[j, 0]:
            visited[j, 0] = True
            queue.append((j, 0))
    while queue:
        j, i = queue.popleft()
        if i == NX - 1:
            return True
        for dj, di in ((-1,0),(1,0),(0,-1),(0,1)):
            nj, ni = j+dj, i+di
            if 0 <= nj < NY and 0 <= ni < NX:
                if not visited[nj, ni] and not mask[nj, ni]:
                    visited[nj, ni] = True
                    queue.append((nj, ni))
    return False

if not bfs_connected(grain_mask):
    sys.exit("ERROR: Pore network is not connected from inlet to outlet. "
             "Adjust grain parameters or seed.")
print("Connectivity check PASSED.")

# ---------------------------------------------------------------------------
# Build row-strip grain regions
# ---------------------------------------------------------------------------
# For each grain cell run, emit one REGION block.
# A "strip" is a contiguous run of grain cells in the same row.
# Region name: grain_r{j:02d}_s{strip_idx:03d}
# Region coordinates span that strip in x, the full cell in y/z.

regions = []   # list of (name, x0, x1, y0, y1)

for j in range(NY):
    in_run  = False
    run_start = 0
    strip_idx = 0
    for i in range(NX):
        if grain_mask[j, i] and not in_run:
            in_run    = True
            run_start = i
        elif not grain_mask[j, i] and in_run:
            in_run = False
            name = f"grain_r{j:02d}_s{strip_idx:03d}"
            x0 = run_start * DX
            x1 = i         * DX
            y0 = j         * DY
            y1 = (j + 1)   * DY
            regions.append((name, x0, x1, y0, y1))
            strip_idx += 1
    if in_run:   # strip runs to last column
        name = f"grain_r{j:02d}_s{strip_idx:03d}"
        x0 = run_start * DX
        x1 = NX        * DX
        y0 = j         * DY
        y1 = (j + 1)   * DY
        regions.append((name, x0, x1, y0, y1))

print(f"Generated {len(regions)} grain region blocks.")

# ---------------------------------------------------------------------------
# Write PFLOTRAN input file
# ---------------------------------------------------------------------------
OUT = "realistic_pore.in"

lines = []
A = lines.append  # shorthand

A("!==============================================================================")
A("! Pore-scale Stokes flow + GIRT tracer transport (pressure-gradient driven).")
A("! Randomly packed circular grains, RSA seed={}, {:d} grains.".format(SEED, len(grains)))
A("! Grid: {:d}x{:d}x{:d}, domain {}x{}x{} m, dx=dy={} m".format(
    NX, NY, NZ, LX, LY, LZ, DX))
A("! Inlet (west): p = {:.1f} Pa;  Outlet (east): p = 0 Pa".format(P_IN))
A("!==============================================================================")
A("")
A("SIMULATION")
A("  SIMULATION_TYPE SUBSURFACE")
A("  PROCESS_MODELS")
A("    SUBSURFACE_FLOW flow")
A("      MODE STOKES")
A("      OPTIONS")
A("        DENSITY   1.d0")
A("        VISCOSITY 1.d0")
A("        ADD_CONVECTION")
A("      /")
A("    /")
A("    SUBSURFACE_TRANSPORT transport")
A("      MODE GIRT")
A("      OPTIONS")
A("        REFERENCE_TEMPERATURE 25.d0")
A("      /")
A("    /")
A("  /")
A("END")
A("")
A("SUBSURFACE")
A("")
A("!------------------------------------------------------------------------------")
A("! Numerical methods")
A("!------------------------------------------------------------------------------")
A("")
A("NUMERICAL_METHODS FLOW")
A("  NEWTON_SOLVER")
A("    RESIDUAL_SCALED_INF_TOL 1.d-2")
A("    ITOL        1.d-8")
A("    ITOL_UPDATE 1.d-2")
A("    RTOL        1.d-8")
A("  /")
A("/")
A("")
A("NUMERICAL_METHODS TRANSPORT")
A("  NEWTON_SOLVER")
A("    ITOL        1.d-8")
A("    ITOL_UPDATE 1.d-2")
A("    RTOL        1.d-8")
A("  /")
A("/")
A("")
A("!------------------------------------------------------------------------------")
A("! Grid: {:d} x {:d} x {:d}, dx = dy = {:.4f} m".format(NX, NY, NZ, DX))
A("!------------------------------------------------------------------------------")
A("")
A("GRID")
A("  TYPE structured")
A("  ORIGIN 0.d0 0.d0 0.d0")
A("  BOUNDS")
A("    0.d0 0.d0 0.d0")
A("    {:.1f}d0 {:.1f}d0 {:.1f}d0".format(LX, LY, LZ))
A("  /")
A("  NXYZ {:d} {:d} {:d}".format(NX, NY, NZ))
A("END")
A("")
A("!------------------------------------------------------------------------------")
A("! Material")
A("!------------------------------------------------------------------------------")
A("")
A("MATERIAL_PROPERTY pore_fluid")
A("  ID 1")
A("  POROSITY 1.0")
A("  TORTUOSITY 1.0")
A("  PERMEABILITY")
A("    PERM_ISO 1.d0")
A("  /")
A("END")
A("")
A("!------------------------------------------------------------------------------")
A("! Fluid / chemistry")
A("!------------------------------------------------------------------------------")
A("")
A("FLUID_PROPERTY")
A("  DIFFUSION_COEFFICIENT 5.d-3")
A("END")
A("")
A("CHEMISTRY")
A("  PRIMARY_SPECIES")
A("    Tracer")
A("  END")
A("  OUTPUT")
A("    TOTAL")
A("    Tracer")
A("  /")
A("END")
A("")
A("!------------------------------------------------------------------------------")
A("! Regions")
A("!------------------------------------------------------------------------------")
A("")
A("REGION all")
A("  COORDINATES")
A("    0.d0 0.d0 0.d0")
A("    {:.1f}d0 {:.1f}d0 {:.1f}d0".format(LX, LY, LZ))
A("  /")
A("END")
A("")
A("! ---- Solid grain row-strip regions ----")
for (name, x0, x1, y0, y1) in regions:
    A(f"REGION {name}")
    A("  COORDINATES")
    A(f"    {x0:.4f}d0 {y0:.4f}d0 0.d0")
    A(f"    {x1:.4f}d0 {y1:.4f}d0 {LZ:.1f}d0")
    A("  /")
    A("END")
    A("")

A("! ---- Boundary face regions ----")
A("REGION inlet_top")
A("  FACE WEST")
A("  COORDINATES")
A("    0.d0 0.8d0 0.d0")
A("    0.d0 1.0d0 0.1d0")
A("  /")
A("END")
A("")
A("REGION inlet_bottom")
A("  FACE WEST")
A("  COORDINATES")
A("    0.d0 0.0d0 0.d0")
A("    0.d0 0.8d0 0.1d0")
A("  /")
A("END")
A("")
A("REGION right")
A("  FACE EAST")
A("  COORDINATES")
A("    {:.1f}d0 0.d0 0.d0".format(LX))
A("    {:.1f}d0 {:.1f}d0 {:.1f}d0".format(LX, LY, LZ))
A("  /")
A("END")
A("")
A("REGION bottom")
A("  FACE SOUTH")
A("  COORDINATES")
A("    0.d0 0.d0 0.d0")
A("    {:.1f}d0 0.d0 {:.1f}d0".format(LX, LZ))
A("  /")
A("END")
A("")
A("REGION top")
A("  FACE NORTH")
A("  COORDINATES")
A("    0.d0 {:.1f}d0 0.d0".format(LY))
A("    {:.1f}d0 {:.1f}d0 {:.1f}d0".format(LX, LY, LZ))
A("  /")
A("END")
A("")
A("!------------------------------------------------------------------------------")
A("! Flow conditions")
A("!------------------------------------------------------------------------------")
A("")
A("FLOW_CONDITION inflow")
A("  TYPE")
A("    LIQUID_PRESSURE DIRICHLET")
A("  /")
A("  LIQUID_PRESSURE {:.1f}d0".format(P_IN))
A("END")
A("")
A("FLOW_CONDITION outflow")
A("  TYPE")
A("    LIQUID_PRESSURE DIRICHLET")
A("  /")
A("  LIQUID_PRESSURE 0.d0")
A("END")
A("")
A("FLOW_CONDITION no_slip")
A("  TYPE")
A("    VELOCITY DIRICHLET")
A("  /")
A("  VELOCITY 0.d0 0.d0 0.d0")
A("END")
A("")
A("!------------------------------------------------------------------------------")
A("! Transport conditions")
A("!------------------------------------------------------------------------------")
A("")
A("TRANSPORT_CONDITION tracer_inlet")
A("  TYPE DIRICHLET_ZERO_GRADIENT")
A("  CONSTRAINT tracer_inlet")
A("    CONCENTRATIONS")
A("      Tracer 1.d0 F")
A("    /")
A("  /")
A("END")
A("")
A("TRANSPORT_CONDITION background")
A("  TYPE DIRICHLET_ZERO_GRADIENT")
A("  CONSTRAINT background")
A("    CONCENTRATIONS")
A("      Tracer 1.d-10 F")
A("    /")
A("  /")
A("END")
A("")
A("!------------------------------------------------------------------------------")
A("! Initial and boundary / prescribed conditions")
A("!------------------------------------------------------------------------------")
A("")
A("INITIAL_CONDITION initial")
A("  FLOW_CONDITION      inflow")
A("  TRANSPORT_CONDITION background")
A("  REGION all")
A("END")
A("")
A("! West: pressure inlet -- top strip injects tracer")
A("PRESCRIBED_CONDITION inlet_top")
A("  FLOW_CONDITION      inflow")
A("  TRANSPORT_CONDITION tracer_inlet")
A("  REGION inlet_top")
A("END")
A("")
A("PRESCRIBED_CONDITION inlet_bottom")
A("  FLOW_CONDITION      inflow")
A("  TRANSPORT_CONDITION background")
A("  REGION inlet_bottom")
A("END")
A("")
A("! East: pressure outlet")
A("PRESCRIBED_CONDITION outlet")
A("  FLOW_CONDITION      outflow")
A("  TRANSPORT_CONDITION background")
A("  REGION right")
A("END")
A("")
A("! No-slip walls")
A("PRESCRIBED_CONDITION bottom_wall")
A("  FLOW_CONDITION      no_slip")
A("  TRANSPORT_CONDITION background")
A("  REGION bottom")
A("END")
A("")
A("PRESCRIBED_CONDITION top_wall")
A("  FLOW_CONDITION      no_slip")
A("  TRANSPORT_CONDITION background")
A("  REGION top")
A("END")
A("")
A("!------------------------------------------------------------------------------")
A("! Material assignment + grain inactivation")
A("!------------------------------------------------------------------------------")
A("")
A("STRATA")
A("  REGION all")
A("  MATERIAL pore_fluid")
A("END")
A("")
for (name, *_) in regions:
    A("STRATA")
    A(f"  REGION {name}")
    A("  INACTIVE")
    A("END")
    A("")

A("!------------------------------------------------------------------------------")
A("! Time stepping")
A("!   Pressure gradient {:.1f} Pa over {:.1f} m => mean u ≈ {:.2f} m/s (open channel)".format(
    P_IN, LX, P_IN / (12. * LX)))
A("!   Pore-network resistance ~2x open channel; u_pore ≈ 0.1 m/s")
A("!   Run to t = 20 s so tracer front fully traverses the domain.")
A("!   Output every 0.2 s => 100 VTK snapshots.")
A("!------------------------------------------------------------------------------")
A("")
A("TIME")
A("  FINAL_TIME            20.d0 s")
A("  INITIAL_TIMESTEP_SIZE 1.d-3 s")
A("  MAXIMUM_TIMESTEP_SIZE 2.d-1 s")
A("END")
A("")
A("!------------------------------------------------------------------------------")
A("! Output")
A("!------------------------------------------------------------------------------")
A("")
A("OUTPUT")
A("  FORMAT VTK")
A("  VARIABLES")
A("    LIQUID_PRESSURE")
A("    STOKES_VELOCITY_X")
A("    STOKES_VELOCITY_Y")
A("    STOKES_VELOCITY_Z")
A("    MATERIAL_ID")
A("  /")
A("  PERIODIC TIME 2.d-1 s")
A("/")
A("")
A("END_SUBSURFACE")

with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"Wrote {OUT} ({len(lines)} lines).")
