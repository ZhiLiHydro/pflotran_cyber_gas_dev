'''
   mapdfn2pflotran.py

   Call methods in mapdfn.py to take output of dfnWorks-Version2.0, create
   equivalent continuous porous medium representation, and write parameters
   (permeability, porosity, tortuosity) to files for use with PFLOTRAN.

   Usage: Edit values for origin, nx, ny, nz, d, k_background, bulk_por,
          tortuosity factor, and h5origin.
          Paths and filenames are hardwired and may also need to be checked.
          As written, they assume script is being called from a subdirectory.
          Then: python mapdfn2pflotran.py

   Dependencies: mapdfn.py
                 numpy
                 h5py

   Author: Emily Stein (ergiamb@sandia.gov)
           Applied Systems Analysis and Research, 8844
           Sandia National Laboratories

   Date: 07/13/18
   SAND Number: SAND2018-7605 O

   Edited by Teresa Portone (tporton@sandia.gov) 11/2020
   to take arguments.

   Edited by Rosie Leone (rleone@sandia.gov) 10/2025 to
   be compatible with dfnWorks-Version2.7


   Assumptions for running:
    1. The center of the DFN is 0,0,0
    2. The files are located in dfn_dir and
      
'''
import os.path as op
import sys, os
sys.path.append(op.abspath(op.dirname(sys.argv[0])))
import numpy as np
from mapdfn import *
from h5py import File
import argparse
from glob import glob


#tolerance to determine if cell size evenly divides domain even if the cell size is not an integer
tolerance = 1e-6
    
def mapdfn2pflotran(DFN,cell_size,matrix_perm,matrix_porosity,tortuosity_factor=0.001,
                    sub_domain = None,sub_origin = None,correction_factor=True,
                    meshing_case="whole_domain"):
    
    print("Starting mapdfn2pflotran!")

    output_dir = "mapdfn_ecpm"
    
    #set up domain first
    if isinstance(cell_size, float) or isinstance(cell_size, int):
       cell_size = cell_size
    elif len(cell_size) != 3:
       error_msg = f"Error: The cell size must be a float or list of three values"
       sys.stderr.write(error_msg)
       sys.exit(1)

    if sub_origin and sub_domain:
      if len(sub_origin) != 3:
         error_msg = f"Error: The sub_origin size must be a list of three values"
         sys.stderr.write(error_msg)
         sys.exit(1)
      if len(sub_domain) != 3:
         error_msg = f"Error: The sub_domain must be a list of three values"
         sys.stderr.write(error_msg)
         sys.exit(1)
      domain = {'x':sub_domain[0],'y':sub_domain[1],'z':sub_domain[2]}
      origin = [-1 * domain['x'] / 2, -1 * domain['y'] / 2, -1 * domain['z'] / 2]
      
      h5origin = sub_origin
      
      filenames = setup_output_dir(output_dir, repo=True)
      
    elif sub_domain or sub_origin:
      error_msg = f"Error: Must specify both sub_domain and sub_origin"
      sys.stderr.write(error_msg)
      sys.exit(1)
      
    else:
      domain = DFN.domain
      # Origin of area to map in DFN domain coordinates (0,0,0 is center of DFN)
      origin = [-1 * domain['x'] / 2, -1 * domain['y'] / 2, -1 * domain['z'] / 2]
      h5origin = [0,0,0]
      filenames = setup_output_dir(output_dir)

    nx, ny, nz, num_cells = setup_domain(domain, cell_size)

    ellipses = readEllipseDFNWorks(DFN.radii,DFN.normal_vectors,DFN.centers)
    

    upscale("./dfnGen_output/",origin,h5origin,ellipses,DFN.aperture,DFN.perm,matrix_perm,matrix_porosity,tortuosity_factor,filenames,nx,ny,nz,cell_size)
    

def main(args):
  parser = argparse.ArgumentParser(description="This script takes the top-level directory of the dfn and maps it to an ecpm, saving the ecpm files in that directory.", 
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("dfn_dir", help="The directory containing the generated dfn files.")
  parser.add_argument("--out_dir", help="The output directory containing the generated ecpm files.", default='mapdfn_ecpm',type=str)
  parser.add_argument( "-d", "--cell_size", help="The cell size to use for the meshing.", default=20, type=float)
  parser.add_argument( "--matrix_perm", help="Permeability of the matrix cells.", default=1e-18, type=float)
  #parser.add_argument( "--effective_diffusivity","--ediff", help="Effective diffusivity.", default=0.000001, type=float)
  parser.add_argument( "--tortuosity_factor", "--tort_factor", help="Tortuosity factor.", default=0.001, type=float)
  parser.add_argument( "--matrix_porosity", help="Porosity of the matrix cells.", default=0.01, type=float)
  parser.add_argument( "--domain", help="Domain length.", nargs=3,type=float)
  parser.add_argument( "--domain_origin", help="Domain origin.", nargs=3,type=float)
  parser.add_argument("--meshing_case", help="A case switch for the meshing to be performed. ",
      choices=['whole_domain', 'repo'], default='whole_domain')

  pa = parser.parse_args(args)

  domain = pa.domain
  domain_origin = pa.domain_origin

  if not op.isdir(pa.dfn_dir): 
    sys.exit("ERROR: you must pass a directory as the dfn_dir argument.")
  
  dfn_dir = op.abspath(pa.dfn_dir)+'/dfnGen_output/'
  if not len(glob(dfn_dir+'radii_Final.dat'))>0:
    sys.exit("ERROR: you must pass a directory holding generated DFN files as the "\
        "dfn_dir argument.")

  """ Setting the properties that are shared no matter the meshing case """
  k_background = pa.matrix_perm
  bulk_por = pa.matrix_porosity
  tortuosity_factor = pa.tortuosity_factor

  """ Setting the properties that depend on the meshing case """
  # origin of area of map in DFN domain coordinates where (0,0,0) is center of DFN
  # PFLOTRAN takes origin at bottom back corner

  domain_origin = [pa.domain_origin[0], pa.domain_origin[1], pa.domain_origin[2]]
  [Lx, Ly, Lz] = [pa.domain[0],pa.domain[1], pa.domain[2]]
  d = pa.cell_size

  if pa.meshing_case=='whole_domain':

    filenames = setup_output_dir(pa.out_dir)
    
    # Origin of area to map in DFN domain coordinates (0,0,0 is center of DFN)
    origin = [-1 * domain[0] / 2, -1 * domain[1] / 2, -1 * domain[2] / 2]
    h5origin = [0,0,0]


  if pa.meshing_case=='repo':

    filenames = setup_output_dir(pa.out_dir,repo=True)
   
    origin = [-1 * domain[0] / 2, -1 * domain[1] / 2, -1 * domain[2] / 2]
    h5origin = domain_origin

  [ nx, ny, nz ] = [ int(Lx/d), int(Ly/d), int(Lz/d) ]

  #if ((Lx % d + Ly % d + Lz % d) > tolerance):
  #      print (Lx%d)
  #      error_msg = f"Error: The cell size you've specified, {d}  m, does not evenly divide the domain. Domain size: {Lx} x {Ly} x {Lz} m^3."
  #      sys.stderr.write(error_msg)
  #      sys.exit(1)
        
  # Call mapdfn functions
  print('Mapping DFN to grid')
  ellipses = readEllipse(dfn_dir+'radii_Final.dat',
      dfn_dir+'normal_vectors.dat',dfn_dir+'translations.dat')
  apertures = readApertures( dfn_dir+'../aperture.dat' )
  perms = readPerms( dfn_dir+'../perm.dat' )

  upscale(dfn_dir,origin,h5origin,ellipses,apertures,perms,k_background,bulk_por,tortuosity_factor,filenames,nx,ny,nz,d)
  
  
  ############################################################

def upscale(dfn_dir,origin,h5origin,ellipses,apertures,perms,k_background,bulk_por,
            tortuosity_factor,filenames,nx,ny,nz,d):

  #ellipses
  fracture = map_dfn(ellipses, origin, nx, ny, nz, d)

  T = apertures * perms

  k_iso = permIso(fracture,T,d,k_background)
  k_aniso = permAniso(fracture,ellipses,T,d,k_background)
  por = porosity(fracture,apertures,d,bulk_por)
  
  ##io
  #arrays for making PFLOTRAN hf file
  x = np.zeros(nx+1,'=f8')
  x[nx] = h5origin[0]+nx*d
  y = np.zeros(ny+1,'=f8')
  y[ny] = h5origin[1]+ny*d
  z = np.zeros(nz+1,'=f8')
  z[nz] = h5origin[2]+nz*d
  a = np.zeros((nx,ny,nz),'=f8')
  khdf5 = np.zeros((nx,ny,nz),'=f8')
  kx = np.zeros((nx,ny,nz),'=f8')
  ky = np.zeros((nx,ny,nz),'=f8')
  kz = np.zeros((nx,ny,nz),'=f8')
  phdf5 = np.zeros((nx,ny,nz),'=f8')

  for k in range(nz):
    z[k] = h5origin[2]+k*d
    for j in range(ny):
      y[j] = h5origin[1]+j*d
      for i in range(nx):
        index = i+nx*j+nx*ny*k
        x[i] = h5origin[0]+i*d
        khdf5[i][j][k] = k_iso[index]
        kx[i][j][k] = k_aniso[index][0]
        ky[i][j][k] = k_aniso[index][1]
        kz[i][j][k] = k_aniso[index][2]
        phdf5[i][j][k] = por[index]
  
  
  # Write isotropic permeability to a gridded dataset for use with PFLOTRAN.
  print ('Writing .h5 file for isotropic permeability field')
  h5file2=File(filenames['isotropic_k'],'w')
  # 3d uniform grid
  h5grp = h5file2.create_group('Permeability')
  # 3D will always be XYZ where as 2D can be XY, XZ, etc. and 1D can be X, Y or Z
  h5grp.attrs['Dimension'] = np.bytes_('XYZ')
  # based on Dimension, specify the uniform grid spacing
  h5grp.attrs['Discretization'] = [d,d,d]
  # again, depends on Dimension
  h5grp.attrs['Origin'] = h5origin
  # leave this line out if not cell centered.  If set to False, it will still
  # be true (issue with HDF5 and Fortran)
  h5grp.attrs['Cell Centered'] = [True]
  h5grp.attrs['Space Interpolation Method'] = np.bytes_('Step')
  h5grp.create_dataset('Data', data=khdf5)
  h5file2.close()
  
  # Write porosity as a gridded dataset for use with PFLOTRAN.
  print ('And also porosity as a gridded dataset')
  h5file2=File(filenames['porosity'],'w')
  # 3d uniform grid
  h5grp = h5file2.create_group('Porosity')
  # 3D will always be XYZ where as 2D can be XY, XZ, etc. and 1D can be X, Y or Z
  h5grp.attrs['Dimension'] = np.bytes_('XYZ')
  # based on Dimension, specify the uniform grid spacing
  h5grp.attrs['Discretization'] = [d,d,d]
  # again, depends on Dimension
  h5grp.attrs['Origin'] = h5origin
  # leave this line out if not cell centered.  If set to False, it will still
  # be true (issue with HDF5 and Fortran)
  h5grp.attrs['Cell Centered'] = [True]
  h5grp.attrs['Space Interpolation Method'] = np.bytes_('Step')
  h5grp.create_dataset('Data', data=phdf5) 
  h5file2.close()
  
  # Write tortuosity as a gridded dataset for use with PFLOTRAN.
  print ('And also tortuosity as a gridded dataset')
  h5file2=File(filenames['tortuosity'],'w')
  # 3d uniform grid
  h5grp = h5file2.create_group('Tortuosity')
  # 3D will always be XYZ where as 2D can be XY, XZ, etc. and 1D can be X, Y or Z
  h5grp.attrs['Dimension'] = np.bytes_('XYZ')
  # based on Dimension, specify the uniform grid spacing
  h5grp.attrs['Discretization'] = [d,d,d]
  # again, depends on Dimension
  h5grp.attrs['Origin'] = h5origin
  # leave this line out if not cell centered.  If set to False, it will still
  # be true (issue with HDF5 and Fortran)
  h5grp.attrs['Cell Centered'] = [True]
  h5grp.attrs['Space Interpolation Method'] = np.bytes_('Step')
  h5grp.create_dataset('Data', data=tortuosity_factor/phdf5) 
  h5file2.close()
  
  # Write anisotropic permeability as a gridded dataset for use with PFLOTRAN.
  print('and anisotropic k')
  h5file3=File(filenames['anisotropic_k'],'w')
  # 3d uniform grid
  h5grp = h5file3.create_group('PermeabilityX')
  # 3D will always be XYZ where as 2D can be XY, XZ, etc. and 1D can be X, Y or Z
  h5grp.attrs['Dimension'] = np.bytes_('XYZ')
  # based on Dimension, specify the uniform grid spacing
  h5grp.attrs['Discretization'] = [d,d,d]
  # again, depends on Dimension
  h5grp.attrs['Origin'] = h5origin
  # leave this line out if not cell centered.  If set to False, it will still
  # be true (issue with HDF5 and Fortran)
  h5grp.attrs['Cell Centered'] = [True]
  h5grp.attrs['Space Interpolation Method'] = np.bytes_('Step')
  h5grp.create_dataset('Data', data=kx) #does this matter that it is also called data?
  
  # 3d uniform grid
  h5grp = h5file3.create_group('PermeabilityY')
  # 3D will always be XYZ where as 2D can be XY, XZ, etc. and 1D can be X, Y or Z
  h5grp.attrs['Dimension'] = np.bytes_('XYZ')
  # based on Dimension, specify the uniform grid spacing
  h5grp.attrs['Discretization'] = [d,d,d]
  # again, depends on Dimension
  h5grp.attrs['Origin'] = h5origin
  # leave this line out if not cell centered.  If set to False, it will still
  # be true (issue with HDF5 and Fortran)
  h5grp.attrs['Cell Centered'] = [True]
  h5grp.attrs['Space Interpolation Method'] = np.bytes_('Step')
  h5grp.create_dataset('Data', data=ky) #does this matter that it is also called data?
  
  # 3d uniform grid
  h5grp = h5file3.create_group('PermeabilityZ')
  # 3D will always be XYZ where as 2D can be XY, XZ, etc. and 1D can be X, Y or Z
  h5grp.attrs['Dimension'] = np.bytes_('XYZ')
  # based on Dimension, specify the uniform grid spacing
  h5grp.attrs['Discretization'] = [d,d,d]
  # again, depends on Dimension
  h5grp.attrs['Origin'] = h5origin
  # leave this line out if not cell centered.  If set to False, it will still
  # be true (issue with HDF5 and Fortran)
  h5grp.attrs['Cell Centered'] = [True]
  h5grp.attrs['Space Interpolation Method'] = np.bytes_('Step')
  h5grp.create_dataset('Data', data=kz) #does this matter that it is also called data?
  h5file3.close()
  
  # Write materials.h5 to inactivate non-fracture cells in PFLOTRAN.
  print('Write material id file for inactivating matrix cells')
  h5file4 = File(filenames['materials'],'w')
  materials_group = h5file4.create_group('Materials')
  iarray = np.zeros((nx*ny*nz),'=i4')
  marray = np.zeros((nx*ny*nz),'=i4')
  for i in range(nx*ny*nz):
    iarray[i] = i+1
    if k_iso[i] == k_background:
      marray[i] = 0
    else:
      marray[i] = 1
  h5dset = materials_group.create_dataset('Cell Ids', data=iarray)
  h5dset = materials_group.create_dataset('Material Ids', data=marray)
  h5file4.close()
  print ('Done!')



def setup_output_dir(output_dir, repo=False):
    """ Create ECPM output directory
    """
    print(f"--> Creating output directory {output_dir}")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    repo = ""
    if repo: repo="rep_"

    filenames = {
        'isotropic_k': output_dir + f'/{repo}isotropic_k.h5',
        'anisotropic_k': output_dir + f'/{repo}anisotropic_k.h5',
        'tortuosity': output_dir + f'//{repo}tortuosity.h5',
        'porosity': output_dir + f'//{repo}porosity.h5',
        'materials': output_dir + f'//{repo}materials.h5'
    }

    return filenames


def setup_domain(domain, cell_size,sub_origin=None):
    """ Initialize domain and discretization

    Parameters
    -----------------
        domain : dictionary
            Domain size in x, y, z from DFN object
        cell_size : list
            Hexahedron cell size in x, y, z

    Returns
    -----------------
        nx : int
            Number of cells in x direction 
        ny : int
            Number of cells in y direction
        nz : int
            Number of cells in z direction 
        num_cells : int 
            Total number of cells in the domain

    Notes
    -----------------
    
    
    
    """
    print("--> Computing discrete domain parameters")
    cell_size = [cell_size, cell_size, cell_size]

    [nx, ny, nz] = [
        int(domain['x'] / cell_size[0]),
        int(domain['y'] / cell_size[1]),
        int(domain['z'] / cell_size[2])
    ]

    num_cells = nx * ny * nz
    

   # if ((domain['x'] % cell_size[0] + domain['y'] % cell_size[1] + domain['z'] % cell_size[2]) > tolerance):
        
   #     error_msg = f"Error: The cell size you've specified, {cell_size} m, does not evenly divide the domain. Domain size: {domain['x']} x {domain['y']} x {domain['z']} m^3."
   #     sys.stderr.write(error_msg)
   #     sys.exit(1)
    print(f"--> Hexahedron edge length {cell_size} m")
    print(f"--> Domain is {nx} x {ny} x {ny} cells. ")
    print(f"--> Total number of cells {num_cells}\n")
    return nx, ny, nz, num_cells




if __name__=='__main__':
  main(sys.argv[1:])
