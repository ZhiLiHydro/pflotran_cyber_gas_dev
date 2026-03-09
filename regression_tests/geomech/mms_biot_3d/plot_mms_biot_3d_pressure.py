import math
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

rho_w = 1000.0
p0 = 101325.0
g_z = -9.80665
z0 = 0.0

workdir = Path(__file__).resolve().parent
h5_candidates = [p for p in workdir.glob("*.h5") if "geomech" not in p.name]
if not h5_candidates:
    raise RuntimeError("No flow HDF5 output found (expected *.h5 without 'geomech').")

h5_path = max(h5_candidates, key=lambda p: p.stat().st_mtime)

with h5py.File(h5_path, "r") as h5:
    time_groups = [k for k in h5.keys() if "Time" in k]
    if not time_groups:
        raise RuntimeError("No time groups found in flow HDF5 output.")
    time_group = sorted(time_groups)[-1]

    group = h5[time_group]
    keys = list(group.keys())

    def find_key(fragment):
        for key in keys:
            if fragment.lower() in key.lower():
                return key
        return None

    pressure_key = find_key("Pressure")
    x_key = find_key("X Coordinate") or find_key("X_Coordinate")
    y_key = find_key("Y Coordinate") or find_key("Y_Coordinate")
    z_key = find_key("Z Coordinate") or find_key("Z_Coordinate")

    if not pressure_key:
        raise RuntimeError("Pressure dataset not found in flow HDF5 output.")
    if not (x_key and y_key and z_key):
        raise RuntimeError("Coordinate datasets not found in flow HDF5 output.")

    pressure = group[pressure_key][...]
    z = group[z_key][...]

p_analytical = p0 + rho_w * g_z * (z - z0)

fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
ax.scatter(p_analytical, pressure, s=18, alpha=0.8)

vmin = min(np.min(p_analytical), np.min(pressure))
vmax = max(np.max(p_analytical), np.max(pressure))
ax.plot([vmin, vmax], [vmin, vmax], "k--", linewidth=1)

ax.set_xlabel("Analytical pressure [Pa]")
ax.set_ylabel("Numerical pressure [Pa]")
ax.set_title("MMS Biot 3D: Pressure (hydrostatic)")

fig.savefig("mms_biot_3d_pressure_compare.png", dpi=150)
