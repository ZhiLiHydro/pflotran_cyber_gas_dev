#!/usr/bin/env python3
"""
animate_realistic.py
====================
Reads all realistic_pore-NNN.vtk snapshots produced by realistic_pore.in
and writes two animated GIFs:

  realistic_pore_flow.gif    — velocity magnitude map with streamlines + particles
  realistic_pore_tracer.gif  — tracer concentration front advancing left → right

Usage
-----
  python animate_realistic.py [output_dir]

  output_dir defaults to the current directory (where the VTK files live).

Dependencies
------------
  numpy, matplotlib, imageio, scipy   (pip install numpy matplotlib imageio scipy)
"""

import sys
import os
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")                     # headless / no display required
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import imageio.v2 as imageio
from scipy.interpolate import RegularGridInterpolator

# ---------------------------------------------------------------------------
# Grid constants — must match realistic_pore.in
# ---------------------------------------------------------------------------
NX, NY   = 80, 40
LX, LY   = 2.0, 1.0
DX, DY   = LX / NX, LY / NY
x_c = np.linspace(DX / 2, LX - DX / 2, NX)   # cell-centre x  (len NX)
y_c = np.linspace(DY / 2, LY - DY / 2, NY)   # cell-centre y  (len NY)

DT_OUTPUT = 0.2   # seconds between consecutive VTK snapshots

# BC_EXCL is computed dynamically below, after reading the first snapshot.
# It is the first column index where the mean fluid-cell ux becomes positive,
# masking the unphysical negative-velocity zone caused by the
# PRESCRIBED_CONDITION pressure traction term at the inlet.
BC_EXCL = None   # set after snapshots are loaded

# ---------------------------------------------------------------------------
# VTK reader
# ---------------------------------------------------------------------------
def read_vtk_scalars(fname):
    """
    Parse a PFLOTRAN ASCII VTK file.

    Returns a dict  {scalar_name: np.ndarray(shape=NX*NY, dtype=float)}.
    When a variable name appears more than once the first occurrence is kept.
    """
    data          = {}
    current_name  = None
    current_vals  = []
    ncells        = NX * NY

    _skip_starts = ("CELL_DATA", "POINT_DATA", "FIELD", "VECTORS",
                    "DATASET", "POINTS", "CELLS", "CELL_TYPES",
                    "DIMENSIONS", "SPACING", "ORIGIN")

    with open(fname) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("SCALARS"):
                # --- flush previous variable ---
                if current_name and current_name not in data:
                    if len(current_vals) >= ncells:
                        data[current_name] = np.array(
                            current_vals[:ncells], dtype=float)
                current_name = line.split()[1]
                current_vals = []

            elif line.startswith("LOOKUP_TABLE"):
                pass                              # skip header row

            elif current_name is not None:
                if not any(line.startswith(k) for k in _skip_starts):
                    current_vals.extend(line.split())

    # flush last variable
    if current_name and current_name not in data:
        if len(current_vals) >= ncells:
            data[current_name] = np.array(current_vals[:ncells], dtype=float)

    return data


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def finish_ax(ax, title, t):
    """Apply uniform axis formatting."""
    ax.set_xlim(0, LX)
    ax.set_ylim(0, LY)
    ax.set_aspect("equal")
    ax.set_xlabel("x  [m]", fontsize=10)
    ax.set_ylabel("y  [m]", fontsize=10)
    ax.set_title(f"{title}     t = {t:.1f} s", fontsize=11, pad=6)
    ax.text(0.01, 0.97, "← inlet (p=4 Pa, tracer=1)", transform=ax.transAxes,
            fontsize=7.5, va="top", color="white",
            bbox=dict(boxstyle="round,pad=0.2", fc="steelblue", alpha=0.7))
    ax.text(0.99, 0.97, "outlet →", transform=ax.transAxes,
            fontsize=7.5, va="top", ha="right", color="white",
            bbox=dict(boxstyle="round,pad=0.2", fc="steelblue", alpha=0.7))


# ---------------------------------------------------------------------------
# Particle tracking helpers
# ---------------------------------------------------------------------------
# Background particles: small, many, colored by initial y
N_INJECT   = 20    # new background particles injected at the inlet per frame
# Path-tracer particles: large yellow, injected in batches
N_PATH              = 12   # particles per batch
PATH_INJECT_FRAMES  = 20   # frames between batch injections (= 2 s at DT_OUTPUT=0.1)

N_SUBSTEPS = 5     # RK4 sub-steps per output interval
BG_SIZE    = 10    # background particle marker size (points²)
PATH_SIZE  = 45    # path-tracer marker size
EXIT_SIZE  = 130   # exit-flash star size

def make_vel_interp(ux_field, uy_field):
    """Build linear interpolators for ux and uy on the cell-centre grid."""
    return (
        RegularGridInterpolator(
            (y_c, x_c), ux_field, method='linear',
            bounds_error=False, fill_value=0.0),
        RegularGridInterpolator(
            (y_c, x_c), uy_field, method='linear',
            bounds_error=False, fill_value=0.0),
    )

def rk4_step(px, py, i_ux, i_uy, dt):
    """Vectorized RK4 step for arrays of particle positions."""
    def vel(x, y):
        pts = np.column_stack([np.clip(y, y_c[0], y_c[-1]),
                               np.clip(x, x_c[0], x_c[-1])])
        return i_ux(pts), i_uy(pts)

    k1x, k1y = vel(px,                   py)
    k2x, k2y = vel(px + 0.5*dt*k1x,     py + 0.5*dt*k1y)
    k3x, k3y = vel(px + 0.5*dt*k2x,     py + 0.5*dt*k2y)
    k4x, k4y = vel(px +    dt*k3x,      py +    dt*k3y)

    return (px + dt*(k1x + 2*k2x + 2*k3x + k4x)/6.0,
            py + dt*(k1y + 2*k2y + 2*k3y + k4y)/6.0)

# Grain mask interpolator — built after reading first VTK snapshot
_grain_interp = None

def build_grain_interp(solid_mask):
    """Build a nearest-neighbour interpolator for the grain (solid) mask."""
    global _grain_interp
    _grain_interp = RegularGridInterpolator(
        (y_c, x_c), solid_mask.astype(float), method='nearest',
        bounds_error=False, fill_value=1.0)   # out-of-bounds → treat as solid

def in_grain(px, py):
    """Boolean array: True where particle (px[i], py[i]) is inside a grain."""
    pts = np.column_stack([np.clip(py, y_c[0], y_c[-1]),
                           np.clip(px, x_c[0], x_c[-1])])
    return _grain_interp(pts) > 0.5

def advance_bg(px, py, pc, i_ux, i_uy):
    """Inject + advance + prune background particles."""
    new_y = np.linspace(DY/2, LY - DY/2, N_INJECT)
    px = np.append(px, np.full(N_INJECT, x_c[BC_EXCL]))
    py = np.append(py, new_y)
    pc = np.append(pc, new_y / LY)

    dt_sub = DT_OUTPUT / N_SUBSTEPS
    for _ in range(N_SUBSTEPS):
        px, py = rk4_step(px, py, i_ux, i_uy, dt_sub)

    keep = ((px >= 0) & (px <= LX) &
            (py >= 0) & (py <= LY) &
            ~in_grain(px, py))
    return px[keep], py[keep], pc[keep]

def advance_path(px, py, pc, i_ux, i_uy):
    """Advance path-tracer particles; return (active, exited_y)."""
    dt_sub = DT_OUTPUT / N_SUBSTEPS
    for _ in range(N_SUBSTEPS):
        px, py = rk4_step(px, py, i_ux, i_uy, dt_sub)

    # Collect y-positions of particles that just crossed the outlet
    exited_y = py[(px > LX) & (py >= 0) & (py <= LY)].copy()

    keep = ((px >= 0) & (px <= LX) &
            (py >= 0) & (py <= LY) &
            ~in_grain(px, py))
    return px[keep], py[keep], pc[keep], exited_y

# Uniform y positions for path-tracer batches
path_inject_y  = np.linspace(DY/2, LY - DY/2, N_PATH)
path_inject_c  = path_inject_y / LY   # color in [0, 1]
# path_inject_x is set after BC_EXCL is computed (below, after snapshots load)


# ---------------------------------------------------------------------------
# Locate VTK files
# ---------------------------------------------------------------------------
workdir  = sys.argv[1] if len(sys.argv) > 1 else "."
pattern  = os.path.join(workdir, "realistic_pore-[0-9]*.vtk")
vtk_files = sorted(glob.glob(pattern))

if not vtk_files:
    sys.exit(
        f"No VTK files found matching:\n  {pattern}\n"
        "Run PFLOTRAN with realistic_pore.in first."
    )

print(f"Found {len(vtk_files)} VTK snapshots — reading …")
snapshots = [read_vtk_scalars(f) for f in vtk_files]
print("  … done.")

# ---------------------------------------------------------------------------
# Pre-compute derived quantities
# ---------------------------------------------------------------------------
mat_id     = snapshots[0]["Material_ID"].reshape(NY, NX)
solid_mask = mat_id < 0.5          # Material_ID==0 => inactive grain cell

# Compute BC_EXCL dynamically: find first column where mean fluid-cell ux > 0.
# The inlet pressure BC (PRESCRIBED_CONDITION) causes unphysical negative
# velocity in a boundary layer several cells wide; mask those columns.
_ux0   = snapshots[-1]["Stokes_Velocity_X"].reshape(NY, NX)
_fluid = ~solid_mask
BC_EXCL = NX // 2   # fallback: half the domain
for _i in range(NX):
    _fl = _fluid[:, _i]
    if _fl.sum() > 0 and _ux0[:, _i][_fl].mean() > 0.0:
        BC_EXCL = _i
        break
print(f"  inlet BC exclusion zone: columns 0..{BC_EXCL-1}  (x < {BC_EXCL*DX:.3f} m)")

# Also mask the inlet BC exclusion zone
solid_mask[:, :BC_EXCL] = True

# Injection x for particles and streamline seeds — first uncontaminated column
path_inject_x = x_c[BC_EXCL]
seed_y   = np.linspace(DY / 2, LY - DY / 2, 60)
seed_pts = np.column_stack([np.full_like(seed_y, path_inject_x), seed_y])

build_grain_interp(solid_mask)     # build once from first snapshot

vmag_all = []
for s in snapshots:
    vm = np.sqrt(s["Stokes_Velocity_X"].reshape(NY, NX) ** 2 +
                 s["Stokes_Velocity_Y"].reshape(NY, NX) ** 2)
    vm[solid_mask] = np.nan
    vmag_all.append(vm)
vmax = float(np.nanpercentile(np.stack(vmag_all), 99))

def snap_vel(s):
    return (np.where(solid_mask, 0.0, s["Stokes_Velocity_X"].reshape(NY, NX)),
            np.where(solid_mask, 0.0, s["Stokes_Velocity_Y"].reshape(NY, NX)))

# Grain overlay image: gray where solid, transparent (NaN) where fluid
grain_overlay = np.where(solid_mask, 0.5, np.nan)   # 0.5 → mid-gray

# ---------------------------------------------------------------------------
# Pre-compute background particle positions (pre-seed domain at t=0)
# ---------------------------------------------------------------------------
print("  pre-seeding background particles …", flush=True)
bg_px, bg_py, bg_pc = np.array([]), np.array([]), np.array([])
ux0, uy0 = snap_vel(snapshots[0])
i_ux0, i_uy0 = make_vel_interp(ux0, uy0)
for _ in range(30):
    bg_px, bg_py, bg_pc = advance_bg(bg_px, bg_py, bg_pc, i_ux0, i_uy0)

all_bg_px = [bg_px.copy()]
all_bg_py = [bg_py.copy()]
all_bg_pc = [bg_pc.copy()]

print(f"  advancing background particles …", flush=True)
for idx in range(1, len(snapshots)):
    ux_s, uy_s = snap_vel(snapshots[idx])
    i_ux, i_uy = make_vel_interp(ux_s, uy_s)
    bg_px, bg_py, bg_pc = advance_bg(bg_px, bg_py, bg_pc, i_ux, i_uy)
    all_bg_px.append(bg_px.copy())
    all_bg_py.append(bg_py.copy())
    all_bg_pc.append(bg_pc.copy())

# ---------------------------------------------------------------------------
# Pre-compute path-tracer particle positions
# ---------------------------------------------------------------------------
print("  advancing path-tracer particles …", flush=True)
pt_px, pt_py, pt_pc = np.array([]), np.array([]), np.array([])

# Inject first batch at frame 0
pt_px = np.append(pt_px, np.full(N_PATH, path_inject_x))
pt_py = np.append(pt_py, path_inject_y)
pt_pc = np.append(pt_pc, path_inject_c)

all_pt_px   = [pt_px.copy()]
all_pt_py   = [pt_py.copy()]
all_pt_pc   = [pt_pc.copy()]
all_exit_y  = [[]]          # exited y-positions per frame

for idx in range(1, len(snapshots)):
    ux_s, uy_s = snap_vel(snapshots[idx])
    i_ux, i_uy = make_vel_interp(ux_s, uy_s)

    if idx % PATH_INJECT_FRAMES == 0:
        pt_px = np.append(pt_px, np.full(N_PATH, path_inject_x))
        pt_py = np.append(pt_py, path_inject_y)
        pt_pc = np.append(pt_pc, path_inject_c)

    pt_px, pt_py, pt_pc, exit_y = advance_path(pt_px, pt_py, pt_pc, i_ux, i_uy)

    all_pt_px.append(pt_px.copy())
    all_pt_py.append(pt_py.copy())
    all_pt_pc.append(pt_pc.copy())
    all_exit_y.append(exit_y.tolist())

print("  … done.")

# ---------------------------------------------------------------------------
# Render frames
# ---------------------------------------------------------------------------
frame_dir = os.path.join(workdir, "anim_frames")
os.makedirs(frame_dir, exist_ok=True)

flow_frames   = []
tracer_frames = []

for idx, (snap, vtk_path) in enumerate(zip(snapshots, vtk_files)):
    t = idx * DT_OUTPUT

    ux = snap["Stokes_Velocity_X"].reshape(NY, NX)
    uy = snap["Stokes_Velocity_Y"].reshape(NY, NX)
    vm = np.where(solid_mask, np.nan, np.sqrt(ux**2 + uy**2))
    tr = np.where(solid_mask, np.nan, snap["Total_Tracer"].reshape(NY, NX))

    ux_s = np.where(solid_mask, 0.0, ux)
    uy_s = np.where(solid_mask, 0.0, uy)

    # ── FLOW frame ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(13, 6.5), dpi=100)
    ax.set_facecolor("#888888")   # gray background → grains show as gray

    im = ax.imshow(vm, origin="lower", extent=[0, LX, 0, LY],
                   cmap="plasma", vmin=0, vmax=vmax,
                   aspect="equal", interpolation="nearest")
    plt.colorbar(im, ax=ax, label="|u|  [m/s]", shrink=0.88, pad=0.02)

    # Streamlines
    try:
        ax.streamplot(
            x_c, y_c, ux_s, uy_s,
            start_points=seed_pts,
            color="white", linewidth=0.6, arrowsize=0.7,
            integration_direction="forward",
            maxlength=8.0,
            zorder=3,
        )
    except Exception as exc:
        print(f"  streamplot warning frame {idx}: {exc}")

    # Background particles — small, cool colormap, colored by initial y
    if len(all_bg_px[idx]) > 0:
        ax.scatter(all_bg_px[idx], all_bg_py[idx],
                   c=all_bg_pc[idx], cmap="cool",
                   vmin=0, vmax=1, s=BG_SIZE,
                   zorder=5, linewidths=0, alpha=0.9)

    # Path-tracer particles — large, colored by initial y
    if len(all_pt_px[idx]) > 0:
        ax.scatter(all_pt_px[idx], all_pt_py[idx],
                   c=all_pt_pc[idx], cmap="autumn",
                   vmin=0, vmax=1, s=PATH_SIZE,
                   zorder=6, linewidths=0.6, edgecolors="white", alpha=0.95)

    # Exit flash — star at outlet for particles that crossed x = LX this frame
    if len(all_exit_y[idx]) > 0:
        ax.scatter(np.full(len(all_exit_y[idx]), LX - DX / 2),
                   all_exit_y[idx],
                   marker="*", s=EXIT_SIZE,
                   c="yellow", edgecolors="orange",
                   linewidths=0.5, zorder=7)

    finish_ax(ax, "Stokes pore-scale flow — random circular grains (velocity magnitude)", t)
    fig.tight_layout()

    fp = os.path.join(frame_dir, f"flow_{idx:03d}.png")
    fig.savefig(fp, bbox_inches="tight", dpi=100)
    plt.close(fig)
    flow_frames.append(imageio.imread(fp))

    # ── TRACER frame ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(13, 6.5), dpi=100)
    ax.set_facecolor("#888888")

    tr_plot = np.where(tr > 0, tr, np.nan)
    im = ax.imshow(tr_plot, origin="lower", extent=[0, LX, 0, LY],
                   cmap="YlOrRd", norm=LogNorm(vmin=1e-10, vmax=1.0),
                   aspect="equal", interpolation="nearest")
    plt.colorbar(im, ax=ax, label="Tracer  [mol/L]", shrink=0.88, pad=0.02)

    finish_ax(ax, "Tracer transport — random circular grain pore network", t)
    fig.tight_layout()

    fp = os.path.join(frame_dir, f"tracer_{idx:03d}.png")
    fig.savefig(fp, bbox_inches="tight", dpi=100)
    plt.close(fig)
    tracer_frames.append(imageio.imread(fp))

    print(f"  frame {idx:3d}/{len(snapshots)-1}  t = {t:.1f} s", flush=True)

# ---------------------------------------------------------------------------
# Assemble GIFs
# ---------------------------------------------------------------------------
def make_durations(n, fps=5, hold_first=1.0, hold_last=2.0):
    d = [1.0 / fps] * n
    if n > 0:
        d[0]  = hold_first
        d[-1] = hold_last
    return d

durs = make_durations(len(flow_frames))

out_flow   = os.path.join(workdir, "realistic_pore_flow.gif")
out_tracer = os.path.join(workdir, "realistic_pore_tracer.gif")

imageio.mimsave(out_flow,   flow_frames,   duration=durs, loop=0)
print(f"Saved  {out_flow}")

imageio.mimsave(out_tracer, tracer_frames, duration=durs, loop=0)
print(f"Saved  {out_tracer}")
