import math
from h5py import *
import numpy as np

prefix = 'distance'

nx = 3
ny = 3
nz = 3
# cell centers for 3x3x3 domain [0,3]^3
x = [0.5, 1.5, 2.5]
nxXny = nx * ny
n = nxXny * nz

# ---- cell-indexed distance (existing 3x3x3 fill test) ----
f = File(prefix + '_cell_indexed.h5', mode='w')
distance = np.zeros(n, dtype=np.double)
for k in range(nz):
    zz = x[k]
    for j in range(ny):
        yy = x[j]
        for i in range(nx):
            xx = x[i]
            distance[i + j * nx + k * nxXny] = math.sqrt(xx * xx + yy * yy + zz * zz)
f.create_dataset('Distance', data=distance)
f.create_dataset('Cell Ids', data=np.arange(1, n + 1, dtype=int))
f.close()

# ---- gridded distance ----
f = File(prefix + '_gridded.h5', mode='w')
h5grp = f.create_group('Distance')
ngx = ngy = ngz = 4
xg = [0., 1., 2., 3.]
h5grp.attrs['Dimension'] = np.bytes_('XYZ')
h5grp.attrs['Discretization'] = [1., 1., 1.]
h5grp.attrs['Origin'] = [0., 0.]
rarray = np.zeros((ngx, ngy, ngz), '=f8')
for k in range(ngz):
    zz = xg[k]
    for j in range(ngy):
        yy = xg[j]
        for i in range(ngx):
            xx = xg[i]
            rarray[i][j][k] = math.sqrt(xx * xx + yy * yy + zz * zz)
h5grp.create_dataset('Data', data=rarray)
f.close()

def porosity_of(xx, yy, zz):
    # valid range over [0,3]^3: 0.10 .. 0.325 at cell centers; 0.10 .. 0.37 at grid nodes
    return 0.10 + 0.02 * xx + 0.03 * yy + 0.04 * zz

# ---- cell-indexed porosity (spatially varying, valid range) ----
f = File('porosity_cell_indexed.h5', mode='w')
porosity = np.zeros(n, dtype=np.double)
for k in range(nz):
    zz = x[k]
    for j in range(ny):
        yy = x[j]
        for i in range(nx):
            xx = x[i]
            porosity[i + j * nx + k * nxXny] = porosity_of(xx, yy, zz)
f.create_dataset('Porosity', data=porosity)
f.create_dataset('Cell Ids', data=np.arange(1, n + 1, dtype=int))
f.close()

# ---- gridded porosity (same formula on node grid as Distance_gridded) ----
f = File('porosity_gridded.h5', mode='w')
h5grp = f.create_group('Porosity')
h5grp.attrs['Dimension'] = np.bytes_('XYZ')
h5grp.attrs['Discretization'] = [1., 1., 1.]
h5grp.attrs['Origin'] = [0., 0.]
parray = np.zeros((ngx, ngy, ngz), '=f8')
for k in range(ngz):
    zz = xg[k]
    for j in range(ngy):
        yy = xg[j]
        for i in range(ngx):
            xx = xg[i]
            parray[i][j][k] = porosity_of(xx, yy, zz)
h5grp.create_dataset('Data', data=parray)
f.close()

def porosity_of_t1(xx, yy, zz):
    # Second time slice: different spatial field, still a valid porosity range.
    # At cell centers: 0.215 .. 0.245; at nodes: 0.20 .. 0.26
    return 0.20 + 0.01 * xx + 0.01 * yy + 0.01 * zz

# ---- gridded porosity time series (Times + multi-slice Data) ----
# Layout matches src/python/gridded_dataset/* : time is the last numpy axis;
# PFLOTRAN Fortran reads dims_h5(1) as time after C/Fortran dim reverse.
nt = 2
times = np.array([0.0, 1.0], dtype=np.float64)
f = File('porosity_gridded_time.h5', mode='w')
h5grp = f.create_group('Porosity')
h5grp.attrs['Dimension'] = np.bytes_('XYZ')
h5grp.attrs['Discretization'] = [1., 1., 1.]
h5grp.attrs['Origin'] = [0., 0.]
h5grp.attrs['Time Units'] = np.bytes_('s')
h5grp.attrs['Time Interpolation Method'] = np.bytes_('LINEAR')
h5grp.attrs['Max Buffer Size'] = [2]
h5grp.create_dataset('Times', data=times)
# shape (ngx, ngy, ngz, nt)
parray_t = np.zeros((ngx, ngy, ngz, nt), '=f8')
for k in range(ngz):
    zz = xg[k]
    for j in range(ngy):
        yy = xg[j]
        for i in range(ngx):
            xx = xg[i]
            parray_t[i][j][k][0] = porosity_of(xx, yy, zz)
            parray_t[i][j][k][1] = porosity_of_t1(xx, yy, zz)
h5grp.create_dataset('Data', data=parray_t)
f.close()

# Expected cell-center values after LINEAR interp at cell centers (node grid).
# DatasetGriddedHDF5InterpolateReal trilinear-interps; for cell centers of the
# unit grid with values defined at nodes, center (0.5,0.5,0.5) averages the
# surrounding nodes — use regression gold from a run rather than hand values
# for gridded centers. Document node extremes for sanity checks.
print('min/max distance cell-indexed:', distance.min(), distance.max())
print('min/max porosity cell-indexed:', porosity.min(), porosity.max())
print('porosity cell 14 (1-based center):', porosity[13])
print('min/max porosity gridded nodes t0:', parray_t[:, :, :, 0].min(),
      parray_t[:, :, :, 0].max())
print('min/max porosity gridded nodes t1:', parray_t[:, :, :, 1].min(),
      parray_t[:, :, :, 1].max())
print('Done')
