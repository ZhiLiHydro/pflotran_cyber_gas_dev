import math
import h5py
import numpy as np
import matplotlib.pyplot as plt

pi = math.pi

h5_path = "mms_elastic_3d-geomech.h5"

with h5py.File(h5_path, "r") as h5:
    time_groups = [k for k in h5.keys() if "Time" in k]
    if not time_groups:
        raise RuntimeError("No time groups found in geomech output.")
    time_group = sorted(time_groups)[-1]

    x = h5["Domain/X"][...]
    y = h5["Domain/Y"][...]
    z = h5["Domain/Z"][...]

    disp_x = h5[f"{time_group}/displacement_x [m]"][...]
    disp_y = h5[f"{time_group}/displacement_y [m]"][...]
    disp_z = h5[f"{time_group}/displacement_z [m]"][...]

# Analytical displacement
f = np.sin(pi * x) * np.sin(pi * y) * np.sin(pi * z)

analytical = {
    "x": f,
    "y": f,
    "z": f,
}

numerical = {
    "x": disp_x,
    "y": disp_y,
    "z": disp_z,
}

fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
for ax, comp in zip(axes, ["x", "y", "z"]):
    ana = analytical[comp]
    num = numerical[comp]
    ax.scatter(ana, num, s=18, alpha=0.8)
    vmin = min(np.min(ana), np.min(num))
    vmax = max(np.max(ana), np.max(num))
    ax.plot([vmin, vmax], [vmin, vmax], "k--", linewidth=1)
    ax.set_xlabel(f"Analytical u_{comp}")
    ax.set_ylabel(f"Numerical u_{comp}")
    ax.set_title(f"Displacement {comp}")

fig.suptitle("MMS Elastic 3D: Numerical vs Analytical", fontsize=12)
fig.savefig("mms_elastic_3d_compare.png", dpi=150)
