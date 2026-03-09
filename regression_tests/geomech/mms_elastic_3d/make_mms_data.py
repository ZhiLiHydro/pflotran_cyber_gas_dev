import math
import h5py
import numpy as np

# Structured grid definition (must match mms_elastic_3d.in)
nx, ny, nz = 4, 4, 4
xmin, ymin, zmin = 0.0, 0.0, 0.0
xmax, ymax, zmax = 1.0, 1.0, 1.0

x = np.linspace(xmin, xmax, nx + 1)
y = np.linspace(ymin, ymax, ny + 1)
z = np.linspace(zmin, zmax, nz + 1)

nnode = (nx + 1) * (ny + 1) * (nz + 1)
cell_ids = np.arange(1, nnode + 1, dtype=np.int32)

# Elastic parameters
E = 1.0e9
nu = 0.25
rho = 1000.0

lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
mu = E / (2.0 * (1.0 + nu))

pi = math.pi

bf_x = np.zeros(nnode)
bf_y = np.zeros(nnode)
bf_z = np.zeros(nnode)

idx = 0
for k, zv in enumerate(z):
    for j, yv in enumerate(y):
        for i, xv in enumerate(x):
            s = math.sin(pi * xv) * math.sin(pi * yv) * math.sin(pi * zv)
            sx = pi * math.cos(pi * xv) * math.sin(pi * yv) * math.sin(pi * zv)
            sy = pi * math.sin(pi * xv) * math.cos(pi * yv) * math.sin(pi * zv)
            sz = pi * math.sin(pi * xv) * math.sin(pi * yv) * math.cos(pi * zv)

            s_xx = -pi * pi * s
            s_yy = -pi * pi * s
            s_zz = -pi * pi * s
            s_xy = (pi * pi) * math.cos(pi * xv) * math.cos(pi * yv) * math.sin(pi * zv)
            s_xz = (pi * pi) * math.cos(pi * xv) * math.sin(pi * yv) * math.cos(pi * zv)
            s_yz = (pi * pi) * math.sin(pi * xv) * math.cos(pi * yv) * math.cos(pi * zv)

            div_x = (lam + 2.0 * mu) * s_xx + (lam + mu) * s_xy + (lam + mu) * s_xz + mu * s_yy + mu * s_zz
            div_y = mu * s_xx + (lam + mu) * s_xy + (lam + 2.0 * mu) * s_yy + (lam + mu) * s_yz + mu * s_zz
            div_z = mu * s_xx + mu * s_yy + (lam + mu) * s_xz + (lam + mu) * s_yz + (lam + 2.0 * mu) * s_zz

            bf_x[idx] = -div_x / rho
            bf_y[idx] = -div_y / rho
            bf_z[idx] = -div_z / rho
            idx += 1

with h5py.File("mms_body_force.h5", "w") as h5:
    h5.create_dataset("Cell Ids", data=cell_ids)
    h5.create_dataset("Body_force_x", data=bf_x)
    h5.create_dataset("Body_force_y", data=bf_y)
    h5.create_dataset("Body_force_z", data=bf_z)


